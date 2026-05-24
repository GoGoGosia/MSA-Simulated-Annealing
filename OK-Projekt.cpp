#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <fstream>

using namespace std;

static const string ALPHABET = "ACGT";
static const int MATCH = 2;
static const int MISMATCH = -1;
static const int GAP_OPEN = -2;
static const int GAP_EXTEND = -1;
static const int GAPGAP = -4; // kara za podwojny gap

struct Alignment {
    vector<vector<int>> M; // Macierz 0 = GAP, 1234 = ACGT
};

// FUNKCJE POMOCNICZE
int charToInt(char c) {
    if (c == 'A') return 1;
    if (c == 'C') return 2;
    if (c == 'G') return 3;
    if (c == 'T') return 4;
    return 0;
}

char intToChar(int x) {
    if (x == 1) return 'A';
    if (x == 2) return 'C';
    if (x == 3) return 'G';
    if (x == 4) return 'T';
    return '-';
}

// PUNKTACJA
int calculateLevenshtein(const string& s1, const string& s2) {
    int n = s1.size();
    int m = s2.size();
    vector<vector<int>> matrix(n + 1, vector<int>(m + 1));

    for (int i = 0; i <= n; i++) matrix[i][0] = i;
    for (int j = 0; j <= m; j++) matrix[0][j] = j;

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int cost;
            if (s1[i - 1] == s2[j - 1]) {
                cost = 0; // Znaki są identyczne - brak kary
            }
            else {
                cost = 1; // Znaki są różne - kara 
            }
            matrix[i][j] = min({ matrix[i - 1][j] + 1, // Usunięcie
                             matrix[i][j - 1] + 1, // Wstawienie
                             matrix[i - 1][j - 1] + cost }); // Zamiana
        }
    }
    return matrix[n][m];
}

string getConsensus(const Alignment& A) {
    string consensus = "";
    int rows = A.M.size();
    int cols = A.M[0].size();

    for (int j = 0; j < cols; j++) {
        int counts[5] = { 0 }; // Licznik -> 0:-, 1:A, 2:C, 3:G, 4:T
        for (int i = 0; i < rows; i++) {
            counts[A.M[i][j]]++;
        }

        int bestNuc = 0;
        int maxCount = -1;

        // Najczęstsza litera w kolumnie
        for (int k = 1; k <= 4; k++) {
            if (counts[k] > maxCount) {
                maxCount = counts[k];
                bestNuc = k;
            }
        }
        if (maxCount > 0) consensus += intToChar(bestNuc);
    }
    return consensus;
}

// GENERATOR INSTANCJI
char randomNuc(mt19937& rng) {
    uniform_int_distribution<int> d(0, 3);
    return ALPHABET[d(rng)];
}

string generateMother(int len, mt19937& rng) {
    string s;
    for (int i = 0; i < len; i++)
        s += randomNuc(rng);
    return s;
}

string mutateWithIndels(const string& m, double p, mt19937& rng) {
    uniform_real_distribution<double> prob(0.0, 1.0);
    string out;

    for (char c : m) {
        double r = prob(rng);

        if (r < p * 0.2) {
            continue; // delecja
        }
        else if (r < p * 0.4) {
            out += randomNuc(rng); // insercja
            out += c;
        }
        else if (r < p) {
            out += randomNuc(rng); // substytucja
        }
        else {
            out += c;
        }
    }
    return out;
}

vector<string> generateInstance(int n, int len, double p) {
    mt19937 rng(time(nullptr));
    vector<string> seqs;

    string mother = generateMother(len, rng);
    seqs.push_back(mother);

    for (int i = 1; i < n; i++)
        seqs.push_back(mutateWithIndels(mother, p, rng));

    return seqs;
}

// ALGORYTMY
Alignment naiveAlign(const vector<string>& seqs) {
    size_t maxLen = 0; // najdłuższa sekwencja
    for (auto& s : seqs)
        maxLen = max(maxLen, s.size());

    Alignment A;
    for (auto& s : seqs) {
        vector<int> row;
        for (char c : s) row.push_back(charToInt(c));
        while (row.size() < maxLen) row.push_back(0); // dopisanie gapów na końcu
        A.M.push_back(row);
    }
    return A;
}

int scoreAlignment(const Alignment& A) {
    int score = 0;
    int rows = A.M.size();
    int cols = A.M[0].size();

    for (int i = 0; i < rows; i++) {
        for (int j = i + 1; j < rows; j++) {

            bool gapInI = false;
            bool gapInJ = false;

            for (int c = 0; c < cols; c++) {

                int a = A.M[i][c];
                int b = A.M[j][c];

                if (a == 0 && b == 0) {
                    score += GAPGAP;
                    gapInI = false;
                    gapInJ = false;
                }
                else if (a == 0) {
                    if (gapInI)
                        score += GAP_EXTEND;
                    else
                        score += GAP_OPEN;
                    gapInI = true;
                    gapInJ = false;
                }
                else if (b == 0) {
                    if (gapInJ)
                        score += GAP_EXTEND;
                    else
                        score += GAP_OPEN;
                    gapInJ = true;
                    gapInI = false;
                }
                else if (a == b) {
                    score += MATCH;
                    gapInI = false;
                    gapInJ = false;
                }
                else {
                    score += MISMATCH;
                    gapInI = false;
                    gapInJ = false;
                }

            }
        }
    }
    return score;
}

Alignment simulatedAnnealing(const Alignment& start, int timeLimit) {
    static mt19937 rng(time(nullptr));
    Alignment best = start;
    Alignment curr = start;

    int bestScore = scoreAlignment(best);
    int currScore = bestScore;

    double Tstart = 20.0;
    double T = Tstart;

    clock_t begin = clock();

    while ((double)(clock() - begin) / CLOCKS_PER_SEC < timeLimit) {
        Alignment cand = curr;
        int rows = cand.M.size();
        int cols = cand.M[0].size();

        int move = rng() % 4; // losowanie sposobu modyfikacji

        // SWAP Zamiana sąsiednich elementów w sekwencji
        if (move == 0) {
            int r = rng() % rows;
            int c = rng() % (cols - 1);
            swap(cand.M[r][c], cand.M[r][c + 1]);
        }

        // Wstawienie przerwy w środek sekwencji i usunięcie ostatniego elementu
        else if (move == 1) {
            int r = rng() % rows;
            int c = rng() % cols;
            if (cand.M[r].back() == 0) {
                cand.M[r].insert(cand.M[r].begin() + c, 0);
                cand.M[r].pop_back();
            }
        }

        // Usunięcie przerwy i dodanie jej na końcu sekwencji
        else if (move == 2) {
            int r = rng() % rows;
            int c = rng() % cols;
            if (cand.M[r][c] == 0) {
                cand.M[r].erase(cand.M[r].begin() + c);
                cand.M[r].push_back(0);
            }
        }

        // Wstawienie przerwy w kolumnie dla części sekwencji 
        else {
            int c = rng() % (cols - 1);
            for (int r = 0; r < rows; r++) {
                if (rng() % 2 && cand.M[r].back() == 0) {
                    cand.M[r].insert(cand.M[r].begin() + c, 0);
                    cand.M[r].pop_back();
                }
            }
        }

        int candScore = scoreAlignment(cand);
        int delta = candScore - currScore;

        if (delta > 0 || exp(delta / T) > (double)rng() / rng.max()) {
            curr = cand;
            currScore = candScore;

            if (currScore > bestScore) {
                best = curr;
                bestScore = currScore;
            }
        }

        double postep = (double)(clock() - begin) / CLOCKS_PER_SEC / timeLimit;
        T = Tstart * (1.0 - postep); // Chłodzenie liniowe: od 20 do 0
    }

    return best;
}

// WYNIKI
void generateFullReport(ostream& out, const Alignment& A, const string& motherDNA,
    int n_seq, int len_dna, double p_mut, int time_sa) {

    int score = scoreAlignment(A);
    string resultConsensus = getConsensus(A);
    int lev = calculateLevenshtein(motherDNA, resultConsensus);

    double similarity = 0;
    if (!motherDNA.empty() && motherDNA.length() > 0) {
        similarity = (1.0 - (double)lev / motherDNA.length()) * 100.0;
    }

    out << "Score: " << score << "\n";
    if (!motherDNA.empty()) {
        out << "Odleglosc Levenshteina:  " << lev << "\n";
        out << "Jakosc rekonstrukcji:    " << similarity << "%\n";
    }
    out << "\n";

    out << "SEKWENCJE:\n";
    if (!motherDNA.empty()) {
        out << "Oryginalne DNA (Mother):\n" << motherDNA << "\n\n";
    }
    out << "Zrekonstruowane DNA (Consensus):\n" << resultConsensus << "\n\n";

    out << "DOPASOWANIE KONCOWE (Alignment):\n";
    for (auto& row : A.M) {
        for (int x : row) out << intToChar(x);
        out << "\n";
    }

    out << "Liczba sekwencji:  " << n_seq << "\n";
    out << "Dlugosc bazowa n:  " << len_dna << "\n";
    if (p_mut < 0) out << "Prawd. mutacji p:  Nieznane (wczytano z pliku)\n";
    else out << "Prawd. mutacji p:  " << p_mut << "\n";
    out << "Czas dzialania SA: " << time_sa << " sek.\n";
}

// OBSUŁUGA PLIKÓW
void saveInstanceToFile(const vector<string>& seqs, const string& filename) {
    ofstream out(filename);
    if (!out) {
        cout << "Blad przy otwieraniu pliku do zapisu!\n";
        return;
    }
    out << seqs.size() << "\n";
    for (const auto& s : seqs)
        out << s << "\n";
    out.close();
    cout << "Instancja zapisana do pliku " << filename << "\n";
}

vector<string> loadInstanceFromFile(const string& filename) {
    ifstream in(filename);
    vector<string> seqs;
    if (!in) {
        cout << "Blad przy otwieraniu pliku do odczytu!\n";
        return seqs;
    }
    int n;
    in >> n;
    in.ignore();
    for (int i = 0; i < n; i++) {
        string s;
        getline(in, s);
        seqs.push_back(s);
    }
    in.close();
    cout << "Instancja wczytana z pliku " << filename << "\n";
    return seqs;
}
void testSA(const Alignment& start, const string& motherDNA) {
    cin.ignore(numeric_limits<streamsize>::max(), '\n');

    long long sumScore = 0;
    int bestScore = -999999;
    int worstScore = 999999;
    int t;
    string temp;

    Alignment bestOverall;

    cout << "Czas SA (sek): ";
    getline(cin, temp);
    if (!temp.empty()) {
        try {
            t = stoi(temp);
        }
        catch (...) {
            cout << "Nieprawidlowy format. Ustawiono domyslna wartosc 10 sek.\n";
        }
    }

    for (int i = 1; i <= 10; i++) {
        Alignment result = simulatedAnnealing(start, t);
        int s = scoreAlignment(result);

        cout << "Proba " << i << ": Score = " << s << endl;

        if (s > bestScore) {
            bestScore = s;
            bestOverall = result; // Zapamiętujemy najlepsze
        }

        if (s < worstScore) {
            worstScore = s;
        }

        sumScore += s;
    }

    double averageScore = (double)sumScore / 10.0;
    string consensus = getConsensus(bestOverall);
    int lev = calculateLevenshtein(motherDNA, consensus);
    double similarity = 0;
    if (!motherDNA.empty()) {
        similarity = (1.0 - (double)lev / motherDNA.length()) * 100.0;
    }

    cout << "\nSREDNI SCORE: " << averageScore << "\n";
    cout << "NAJLEPSZY SCORE: " << bestScore << "\n";
    cout << "NAJGORSZY SCORE: " << worstScore << "\n";
    cout << "CONSENSUS: " << consensus << "\n";
    cout << "\nLEVENSHTEIN: " << lev << "\n";
    cout << "REKONSTRUKCJA: " << similarity << "%\n";
}


int main() {
    vector<string> seqs;
    string motherDNA = "";
    Alignment naive, meta;

    int n = 4;
    int len = 400;
    double p = 0.05;
    int t = 10;

    while (true) {
        cout << "\n--- MULTIPLE SEQUENCE ALIGNMENT ---\n";
        cout << "1. Generuj instancje\n";
        cout << "2. Algorytm naiwny\n";
        cout << "3. Metaheurystyka (SA)\n";
        cout << "4. Wyswietl wynik\n";
        cout << "5. Zapisz instancje do pliku\n";
        cout << "6. Wczytaj instancje z pliku\n";
        cout << "7. Zapisz wynik do pliku\n";
        cout << "8. Testy metaheurystyki (10 prob)\n";
        cout << "0. Wyjscie\n\n> ";

        int op;
        cin >> op;

        switch (op) {
        case 0:
            return 0;

        case 1: {
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            string temp;

            cout << "\nLiczba sekwencji: ";
            getline(cin, temp);
            if (!temp.empty()) {
                try {
                    n = stoi(temp);
                }
                catch (...) {
                    cout << "Nieprawidlowy format. Ustawiono domyslna wartosc 4.\n";
                }
            }

            cout << "Dlugosc sekwencji: ";
            getline(cin, temp);
            if (!temp.empty()) {
                try {
                    len = stoi(temp);
                }
                catch (...) {
                    cout << "Nieprawidlowy format. Ustawiono domyslna wartosc 400.\n";
                }
            }

            cout << "P(mutacji): ";
            getline(cin, temp);
            if (!temp.empty()) {
                try {
                    replace(temp.begin(), temp.end(), ',', '.');
                    double val = stod(temp);
                    if (val < 0.0 || val > 1.0) {
                        cout << "Nieprawidłowa wartość. Ustawiono domyslna wartosc 0.05.\n";

                    }
                    else {
                        p = val;
                    }
                }
                catch (...) {
                    cout << "Nieprawidlowy format. Ustawiono domyslna wartosc 400.\n";
                }
            }

            seqs = generateInstance(n, len, p);
            if (!seqs.empty()) motherDNA = seqs[0];
            cout << "Instancja wygenerowana.\n";
            break;
        }

        case 2:
            naive = naiveAlign(seqs);
            cout << "Score naiwny: " << scoreAlignment(naive) << "\n";
            break;

        case 3: {
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            string temp;

            cout << "Czas SA (sek): ";
            getline(cin, temp);
            if (!temp.empty()) {
                try {
                    t = stoi(temp);
                }
                catch (...) {
                    cout << "Nieprawidlowy format. Ustawiono domyslna wartosc 10 sek.\n";
                }
            }

            meta = simulatedAnnealing(naive, t);
            cout << "Score metaheurystyki: " << scoreAlignment(meta) << "\n";
            break;
        }

        case 4:
            generateFullReport(cout, meta, motherDNA, n, len, p, t);
            break;

        case 5: {
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            string fname;
            cout << "Podaj nazwe pliku: ";
            getline(cin, fname);
            if (fname.empty()) {
                fname = "instancja.txt";
            }
            saveInstanceToFile(seqs, fname);
            break;
        }

        case 6: {
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            string fname;
            cout << "Podaj nazwe pliku: ";
            getline(cin, fname);
            if (fname.empty()) {
                cout << "Nie podano nazwy pliku!\n";
            }
            else {
                seqs = loadInstanceFromFile(fname);
                if (!seqs.empty()) motherDNA = seqs[0];
            }
            break;
        }

        case 7: {
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            string fname;
            cout << "Podaj nazwe pliku: ";
            getline(cin, fname);
            if (fname.empty()) {
                fname = "wynik.txt";
            }
            ofstream out(fname);
            if (out) {
                generateFullReport(out, meta, motherDNA, n, len, p, t);
                out.close();
                cout << "Wynik zapisany do pliku " << fname << "\n";
            }
            break;
        }
        case 8: {
            if (seqs.empty()) {
                cout << "Najpierw wygeneruj lub wczytaj instancje!\n";
            }
            else {
                naive = naiveAlign(seqs);
                testSA(naive, motherDNA);
            }
            break;
        }

        default:
            cout << "Nieprawidlowy wybor!\n";
            break;
        }
    }

    return 0;
}