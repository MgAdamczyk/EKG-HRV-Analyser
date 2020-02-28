#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <libalglib/interpolation.h>
#include <libalglib/fasttransforms.h>
#include "libalglib/stdafx.h"
#include <stdlib.h>
#include <chrono>

using namespace std::chrono;
using namespace std;
using namespace alglib;

bool getFileContent(string fileName, vector<float> & vecOfStrs)
{

    // Open the File
    ifstream in(fileName.c_str());

    // Check if object is valid
    if(!in)
    {
        cerr << "Cannot open the File : " << fileName << endl;
        return false;
    }

    string str;
    // Read the next line from File untill it reaches the end.
    while (getline(in, str))
    {
        // Line contains string of length > 0 then save it in vector
        if(str.size() > 0)
            vecOfStrs.push_back(stof( str ));
    }
    // Close The File
    in.close();

    return true;
}


float getSignalLength(const vector<float> &signal)
{
    /// Obliczanie czasu trwania sygnału ///
    float  endVec = 0;
    for (int i = 0; i<signal.size(); i++)
    {
        endVec += signal[i];
    }
    return endVec;
}

vector<float> convertToMs(const vector<float> &data, const int N)
{
    /// Przeliczanie danych na milisekundy ///
    vector<float> dataInMs;

    for (int i = 0; i < N; i++) {
        dataInMs.push_back(data[i] * 1000);
    }
    return dataInMs;
}

vector<float> accumulateData(const vector<float> &data, const int N)
{
    /// Obliczanie skumulowanych danych w sekundach [s] ///
    vector<float> dataAccumulated;
    float sum = 0;

    for( int i = 0; i < N; i++ )
    {
        sum += data[i];
        dataAccumulated.push_back(sum);
    }

    return dataAccumulated;
}

vector<float> cubicSpline(const vector<float> x, const vector<float> y, int nframe)
{
    alglib::spline1dinterpolant s;
    alglib::ae_int_t natural_bound_type = 2;

    double *newx = new double[nframe];
    double *newy = new double[nframe];

    for( int i = 0; i < nframe; i++)
    {
        newx[i] = x[i];
        newy[i] = y[i];
    }

    alglib::real_1d_array *rx = new alglib::real_1d_array();
    rx->setcontent(nframe, newx);
    alglib::real_1d_array *ry = new alglib::real_1d_array();
    ry->setcontent(nframe, newy);

    alglib::spline1dbuildcubic(
            *const_cast<const alglib::real_1d_array*>(rx),
            *const_cast<const alglib::real_1d_array*>(ry),
            s);

    float* sigI = new float[nframe];

    for(int i = 0; i<nframe; i++)
    {
        sigI[i] = alglib::spline1dcalc(s, i*(getSignalLength(y) / nframe));
    }

    vector<float> sigIVector;

    for(int i = 0; i < nframe; i++)
    {
        sigIVector.push_back(sigI[i]);
    }

    return sigIVector;
}

float calculateMeanRR(const vector<float> &loadedVecRPeaks, const int N)
{
    /// Obliczanie meanRR ///
    float meanRR = 0;

    for(int i = 0; i < N; i++)
    {
        meanRR += loadedVecRPeaks[i];
    }

    return meanRR / N;
}

float calculateSdn(const vector<float> &loadedVecRPeaks, const int N, float &meanRR)
{
    /// Obliczanie SDNN ///
    float sdn = 0;

    for (int i = 0; i < N; i++)
    {
        sdn = sdn + pow(loadedVecRPeaks[i] - meanRR, 2);
    }
    sdn = sqrt(sdn/N-1);

    return sdn;
}

float calculateRmssd(const vector<float> &loadedVecRPeaks, const int N)
{
    /// Obliczanie RMSSD ///
    float rmssd = 0;
    float sumSquare = 0;

    for (int i = 0; i < N-1; i++)
    {
        sumSquare = sumSquare + pow(loadedVecRPeaks[i+1] - loadedVecRPeaks[i],2);
    }
    rmssd = sqrt(sumSquare / (N-1));

    return rmssd;
}

float calculateSda(const vector<float> &loadedVecRPeaks, const int nRange, float *average)
{
    /// Obliczanie SDANN ///
    float sda = 0;
    float suma = 0;
    int nrPrzedzialu = 1;
    float bufor = 0; //tj. odliczanie 5-min czasu
    int liczba = 0;
    int i = 0;

    while(nrPrzedzialu < nRange+1)
    {
        if(bufor > (5*60*1000))
        {
            float avr = suma / liczba;
            average[nrPrzedzialu-1] = avr;
            liczba = 0;
            suma = 0;
            bufor = bufor - 5*60*1000;
            nrPrzedzialu = nrPrzedzialu + 1;
        }
        suma += loadedVecRPeaks[i];
        bufor += loadedVecRPeaks[i];
        liczba++;
        i++;
    }

    // średnia ze wszystkich średnich w przedziałach
    float sum_avr = 0;

    for (int i = 0; i < nRange; i++)
    {
        sum_avr = sum_avr + average[i];
    }
    float mainAvr = sum_avr / nRange;

    // odchylenie standardowe średnich w przedziałach od średniej mainAvr
    for (int i = 0; i < nRange-1; i++)
    {
        sda += pow(average[i] - mainAvr,2);
    }
    sda = sqrt(sda/(nRange-1));

    return sda;
}

float calculateSdi(const vector<float> &loadedVecRPeaks, const int nRange, float *average)
{
    /// Obliczanie SDNN index ///
    float sdi = 0;
    float deviation[nRange];

    float suma = 0;
    int nrPrzedzialu = 1;
    float bufor = 0; //tj. odliczanie 5 min czasu
    int liczba = 0;
    int i = 0;

    while (nrPrzedzialu < nRange+1)
    {
        if(bufor > 5*60*1000)
        {
            float dev = sqrt(suma/(liczba-1));
            deviation[nrPrzedzialu-1] = dev;
            liczba = 0;
            suma = 0;
            bufor = bufor - 5*60*1000;
            nrPrzedzialu++;
            if(nrPrzedzialu == nRange + 1) break;
        }
        suma += pow(loadedVecRPeaks[i] - average[nrPrzedzialu-1],2);
        bufor += loadedVecRPeaks[i];
        liczba++;
        i++;
    }

    // śrdenia ze wszystkich odchyleń uzskanych  w przedziałach
    suma = 0;
    for (int i = 0; i < nRange; i++) suma += deviation[i];

    sdi = suma / nRange;

    return sdi;
}

float calculateNn50(const vector<float> &loadedVecRPeaks, const int N)
{
    /// Obliczanie nn50 ///
    float nn50 = 0;
    for (int i = 0; i < N-1; i++)
    {
        float diff = abs(loadedVecRPeaks[i+1] - loadedVecRPeaks[i]);
        if(diff > 50) nn50++;
    }

    return nn50;
}

float calculateP50nn(const float nn50, const int N)
{
    /// Obliczanie p50nn ///
    return nn50 / N * 100;
}


vector<float> calculatePowerSpectrum(const vector<float> &newValuesY, const int N)
{
    /// Obliczanie FFT ///
    double *input = new double[N];

    for(int i = 0; i < N; i++)
    {
        input[i] = newValuesY[i];
    }

    real_1d_array y;
    y.setcontent(N, input);
    complex_1d_array fftRes;
    fftr1d(y, fftRes);

    /// Obliczanie mocy widma ///
    vector<float> powerSpectrum;

    for(int i = 0; i < N/2; i++)
    {
        if(i == 0) powerSpectrum.push_back(pow(sqrt(fftRes[i].x * fftRes[i].x + fftRes[i].y * fftRes[i].y) , 2 ) / N);
        else powerSpectrum.push_back(2 * pow(sqrt(fftRes[i].x * fftRes[i].x + fftRes[i].y * fftRes[i].y) , 2 ) / N);
    }

    return powerSpectrum;
}

vector<float> calculateFreqSpectrum(const float step, const int N)
{
    /// Obliczanie wektora częstotliwości ///
    vector<float> freqSpectrum;
    float fp = 1/step;
    float step_f = fp/N;

    for(int i = 0; i < N/2; i++)
    {
        freqSpectrum.push_back(i*step_f);
    }

    return freqSpectrum;
}


/// MAIN ///

int main() {

    auto start = high_resolution_clock::now();

    /// Ścieżka do pliku z danymi testowymi ///
    string fileName = "/home/gosia/ImplFunkcjonalnosci/nsr001.dat";

    /// Inicjalizacja parametrów wynikowych ///
    // Parametry statystyczne
    float meanRR = 0;
    float sdn = 0;
    float sda = 0;
    float rmssd = 0;
    float sdi = 0;
    float nn50 = 0;
    float p50nn = 0;

    //Parametry częstotliwościowe
    float TP = 0;
    float HF = 0;
    float LF = 0;
    float VLF = 0;
    float ULF = 0;
    float LFHF = 0;

    //Do wykresu
    vector<float> freqSpectrum;
    vector<float> powerSpectrum;

    /// Wektor z danymi testowymi z pliku ///
    vector<float> loadedVec; //wczytane dane w [s]
    vector<float> loadedVecRPeaks; //potem przeliczone na [ms] z loadedVec

    /// Wczytanie daych z pliku do testów ///
    // informacja o rezultacie wczytywania danych z pliku
    bool result = getFileContent(fileName, loadedVec);

    if (result) cout << "Udało się wczytać plik." << endl << endl;
    else {
        cout << "Nie udało się wczytać pliku!!!" << endl << endl;
        return -1;
    }

    /// TEST: Wypisanie wektora odczytanego z pliku ///
    /*
    for( int i = 0; i < loadedVec.size(); i++ )
    {
        cout << loadedVec[i] << endl;
    }
    */

    /// Wyznaczenie liczby pików R w sygnale ///
    int N = loadedVec.size();

    loadedVecRPeaks = convertToMs(loadedVec, N);
    vector<float> dataAccumulated = accumulateData(loadedVec, N);

    meanRR = calculateMeanRR(loadedVecRPeaks, N);
    sdn = calculateSdn(loadedVecRPeaks, N, meanRR);
    rmssd = calculateRmssd(loadedVecRPeaks, N);

    /// Obliczanie ilości przedziałów pięciominutowych w sygnale ///
    float a = getSignalLength(loadedVecRPeaks);
    int nRange = floor(a / (5 * 60 * 1000));

    float average[nRange];

    sda = calculateSda(loadedVecRPeaks, nRange, average);
    sdi = calculateSdi(loadedVecRPeaks, nRange, average);
    nn50 = calculateNn50(loadedVecRPeaks, N);
    p50nn = calculateP50nn(nn50, N);

    /// Interpolacja tachogramu ///
    vector<float> newValuesY = cubicSpline(dataAccumulated, loadedVec, N);

    cout << "Obliczanie wektora mocy widma... ";
    powerSpectrum = calculatePowerSpectrum(newValuesY, N);
    cout << "Gotowe." << endl << endl;

    /// Obliczanie kroku czasowego w wyniku interpolacji ///
    float step = getSignalLength(loadedVec) / N;

    cout << "Obliczanie wektora częstotliwości... ";
    freqSpectrum = calculateFreqSpectrum(step, N);
    cout << "Gotowe." << endl << endl;

    /// ############################## ///
    /// Zapis wektorów do pliku ///
    cout << "Zapisywanie wektorów do pliku... ";

    ofstream zapis("/home/gosia/ImplFunkcjonalnosci/wyniki.txt");

    zapis << "Wektor częstotliwości.\n";

    for(int i = 0; i < N/2; i++)
    {
        zapis << freqSpectrum[i] << ";";
    }

    zapis << "\n\nWektor mocy widma.\n";

    for(int i = 0; i < N/2; i++)
    {
        zapis << powerSpectrum[i] << ";";
    }

    zapis.close();
    cout << "Gotowe." << endl << endl;
    /// ############################## ///

    /// Obliczanie paramterów częstotliwościowych ///
    float fp = 1/step;
    float step_f = fp/N;

    for(int i = 0; i< N/2; i++)
    {
        if(0.4>=i*step_f && i*step_f>0.15)
        {
            HF += powerSpectrum[i];
        }
        else if(0.15>=i*step_f && i*step_f>0.04)
        {
            LF += powerSpectrum[i];
        }
        else if(0.04>=i*step_f && i*step_f>0.003)
        {
            VLF += powerSpectrum[i];
        }
        else if(0.003>=i*step_f)
        {
            ULF += powerSpectrum[i];
        }
    }
    LFHF = LF/HF;
    TP = HF + LF + VLF + ULF;

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Czas wykonywania programu: " << duration.count() << " us" << endl << endl;

    /// Wypisywanie obliczonych parametrów ///

    cout << "meanRR : " << meanRR << endl;
    cout << "sdn : " << sdn << endl;
    cout << "sda : " << sda << endl;
    cout << "rmssd : " << rmssd << endl;
    cout << "sdi : " << sdi << endl;
    cout << "nn50 : " << nn50 << endl;
    cout << "p50nn : " << p50nn << endl;


    cout<<endl<<endl;

    cout << "TP : " << TP << endl;
    cout << "HF : " << HF << endl;
    cout << "LF : " << LF << endl;
    cout << "VLF : " << VLF << endl;
    cout << "ULF : " << ULF << endl;
    cout << "LFHF : " << LFHF << endl;

    return 0;
}
