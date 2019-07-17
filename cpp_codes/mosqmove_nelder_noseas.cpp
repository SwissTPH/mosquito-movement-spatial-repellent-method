# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <fstream>

# include "asa047.hpp"


using namespace std;

int main ( );
void test01 ( );
double mosqmove ( double x[32] );


//****************************************************************************

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for mosqmove
//       - uses the asa047 library
//

{
  timestamp ( );
  cout << "\n";
 
  test01 ( );

//
//  Terminate.
//
  cout << "\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************

void test01 ( )

{

// for nelmin
  int i;
  int icount;
  int ifault;
  int kcount;
  int konvge;
  int n;
  int numres;
  double reqmin;
  double *start;
  double *step;
  double *xmin;
  double ynewlo;

  n = 32;

  start = new double[n];
  step = new double[n];
  xmin = new double[n];

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Apply NELMIN to mosqmove function.\n";

   start[0] = -1;
   start[1] = -1;
   start[2] = -1;
   start[3] = -1;
   start[4] = -1;
   start[5] = -1;
   start[6] = -1;
   start[7] = -1;
   start[8] = -1;
   start[9] = -1;
   start[10] = -1;
   start[11] = -1;
   start[12] = -1;
   start[13] = -1;
   start[15] = -1;
   start[16] = -1;
   start[17] = -1;
   start[18] = -1;
   start[19] = -1;
   start[20] = -1;
   start[21] = -1;
   start[22] = -1;
   start[23] = -1;
   start[24] = -1;
   start[25] = -1;
   start[26] = -1;
   start[27] = -1;
   start[28] = -1;
   start[29] = -2;
   start[30] = -0.01;
   start[31] = -0.01;

  reqmin = 1.0E-08;

  step[0] = 1.0;
  step[1] = 1.0;
  step[2] = 1.0;
  step[3] = 1.0;
  step[4] = 1.0;
  step[5] = 1.0;
  step[6] = 1.0;
  step[7] = 1.0;
  step[8] = 1.0;
  step[9] = 1.0;
  step[10] = 1.0;
  step[11] = 1.0;
  step[12] = 1.0;
  step[13] = 1.0;
  step[14] = 1.0;
  step[15] = 1.0;
  step[16] = 1.0;
  step[17] = 1.0;
  step[18] = 1.0;
  step[19] = 1.0;
  step[20] = 1.0;
  step[21] = 1.0;
  step[22] = 1.0;
  step[23] = 1.0;
  step[24] = 1.0;
  step[25] = 1.0;
  step[26] = 1.0;
  step[27] = 1.0;
  step[28] = 1.0;
  step[29] = 1.0;
  step[30] = 1.0;
  step[31] = 1.0;

  
  konvge = 10;
  kcount = 500;

  cout << "\n";
  cout << "  Starting point X:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << start[i] << "\n";
  }

  ynewlo = mosqmove ( start );

  cout << "\n";
  cout << "  F(X) = " << ynewlo << "\n";

  nelmin ( mosqmove, n, start, xmin, &ynewlo, reqmin, step,
    konvge, kcount, &icount, &numres, &ifault );

  cout << "\n";
  cout << "  Return code IFAULT = " << ifault << "\n";
  cout << "\n";
  cout << "  Estimate of minimizing value X*:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << xmin[i] << "\n";
  }

  cout << "\n";
  cout << "  F(X*) = " << ynewlo << "\n";

  cout << "\n";
  cout << "  Number of iterations = " << icount << "\n";
  cout << "  Number of restarts =   " << numres << "\n";

  delete [] start;
  delete [] step;
  delete [] xmin;

  return;
}



// This method uses the nelder-mead optimization procedure

double mosqmove ( double x[32] )

{

// initialise variables 
double loglik = 0.0;
// for village 1 obs=30*72=2160 (one row per house per day)
int house[2160]; // house number
float east[2160]; // house coordinate - longitudes
double south[2160]; // house coordinates - latitudes
int exptDay[2160]; // experimental day
int obsmosq[2160]; // number of observed/collected mosquitoes
int combo[2160]; // whether the household has a coil (combo=1) or not (combo=0)
int notmissing[2160];

int numHouses=30; // total number of houses in the village

double pr[30]; // baseline distribution of mosquitoes per house

// input parameters
// input baseline distributions of mosquitoes in each house (pr)
pr[0] = exp(x[0]); pr[1] = exp(x[1]); pr[2] = exp(x[2]); pr[3] = exp(x[3]); pr[4] = exp(x[4]); pr[5] = exp(x[5]);
pr[6] = exp(x[6]); pr[7] = exp(x[7]); pr[8] = exp(x[8]); pr[9] = exp(x[9]); pr[10] = exp(x[10]); pr[11] = exp(x[11]);
pr[12] = exp(x[12]); pr[13] = exp(x[13]); pr[14] = exp(x[14]); pr[15] = exp(x[15]); pr[16] = exp(x[16]); pr[17] = exp(x[17]);
pr[18] = exp(x[18]); pr[19] = exp(x[19]); pr[20] = exp(x[20]); pr[21] = exp(x[21]); pr[22] = exp(x[22]); pr[23] = exp(x[23]);
pr[24] = exp(x[24]); pr[25] = exp(x[25]); pr[26] = exp(x[26]); pr[27] = exp(x[27]); pr[28] = exp(x[28]); pr[29] = 0.3;


double lambda = exp(x[29]); // mean distance moved between houses by the repelled mosquitoes
double spatialRepFactor = exp(x[30]) / (1 + exp(x[30]));  // the proportion of mosquitoes diverted from households using repellents
double hseFactor = exp(x[31]) / (1 + exp(x[31])); // the proportion of mosquitoes diverted elsewhere as opposed to households

// scale prs to one and make a vector of the baseline proportions (prs scaled to one)
double sump = 0;
for (int h=0; h<numHouses; h++) {
     sump = sump + pr[h];
}
for (int h=0; h<numHouses; h++) {
   pr[h] = pr[h]/sump;  
}

// read in data for which to estimate the parameters
    double x1;
    ifstream inFile1;
    inFile1.open("aggTrapSimul.txt");  // dataset as a text file

	if (!inFile1) {
        cout << "Unable to open file";
        exit(1); // terminate with error
   }
   	int i1 = 0;
	int j1 = 0;
	int k1 = 0;
    while (inFile1 >> x1) {
        k1 = i1%7;
        j1 = int(i1/7);
		if (k1==0) (house[j1] = (x1-1));
		if (k1==1) (east[j1] = x1);
        if (k1==2) (south[j1] = x1);
		if (k1==3) (exptDay[j1] = x1);
        if (k1==4) (obsmosq[j1] = x1);
		if (k1==5) (combo[j1] = x1);
        if (k1==6) (notmissing[j1] = x1);
        i1++;
	}
	inFile1.close();

// set up every pair of houses in the villages

int house1[30*30];
int house2[30*30];
double baseph1[30*30];
double south1[30*30];
double south2[30*30];
double east1[30*30];
double east2[30*30];
int combo1[30*30];
int combo2[30*30];
double distance[30*30]; // distance between any two households
double outgoingh1[30*30]; // mosquitoes outgoing from house h1 using spatial repellent
double probGoToh2[30*30]; // probability of moving from house h1 to house h2
double sumProbGoToh2[30]= {0.0}; // sum of probabilities of mosquitoes moving to house h2
double incomingh2[30*30]={0.0}; // mosquitoes incoming to house h2
double incomingh[30] = {0.0}; // total mosquitoes incoming to house h
double outgoingh[30]; // total mosquitoes outgoing from house h
double predph[30]; // predicted mosquitoes in house h
double dissaph[30]; // total mosquitoes from house h diverted elsewhere
double adjdissaph[30]; // adjusted total mosquitoes from house h diverted elsewhere
double adjpredph[30]; // adjusted predicted mosquitoes in house h
double sumpredph = 0.0;  // sum of predicted mosquitoes per household
double sumincoming = 0.0; // sum of all incomning mosquitoes
double sumoutgoing = 0.0; // sum of all outgoing mosquitoes
double sumdiss = 0.0; // sum of all mosquitoes diverted elsewhere
double ntotal = 1000; // total mosquitoes that exist from all households
double sumobsmosq = 0.0; // sum of all observed mosquitoes
double dissapear; // difference between total mosquitoes that exist and those observed
double obsmosqa[31]; // simulated total mosquitoes including the mosquitoes diverted elsewhere
double adjpredpha[31]; // adjusted predicted mosquitoes in house h accounting for mosquitoes diverted elsewhere

int m;
int numDays=72;

// FOR EACH EXPERIMENTAL DAY
for (int d=0; d<numDays; d++) {
  
 // sets up house h1 for the pairs
 for (int h1=0; h1<(30*30); h1++) {
    m = (h1/30);
	house1[h1] = m;
	baseph1[h1] = pr[m];
	south1[h1] = south[m];
	east1[h1] = east[m];
	combo1[h1] = combo[m + (d*30)];
 }


 // sets up house h2
 for (int h2=0; h2<(30*30); h2++) {
    m = (h2 % 30);
	house2[h2] = m;
	south2[h2] = south[m];
	east2[h2] = east[m];
	combo2[h2] = combo[m + (d*30)];
	}

 // calculate distance and probability of movement
for (int h1=0; h1<(30*30); h1++) {
	// calculate harvesian distances between households
    distance[h1] = 6367 * (2 * atan2(sqrt(pow(sin(((east2[h1] - east1[h1]) * 0.0174532925199433)/2.0), 2) + cos(east1[h1]*0.0174532925199433) * cos(east2[h1]*0.0174532925199433) * pow(sin(((south2[h1] - south1[h1]) * 0.0174532925199433)/2.0), 2)), sqrt(1-(pow(sin(((east2[h1] - east1[h1]) * 0.0174532925199433)/2.0), 2) + cos(east1[h1]*0.0174532925199433) * cos(east2[h1]*0.0174532925199433) * pow(sin(((south2[h1] - south1[h1]) * 0.0174532925199433)/2.0), 2)))));

    // outgoing mosquitoes from h1
	outgoingh1[h1] = baseph1[h1] * (spatialRepFactor*static_cast<double>(combo1[h1]));
	
 	// assume moving mosquitoes do not end up in house with repellant
	// the probabiity of moving from h1 to h2 is a function of the distance between h1 and h2 and presence of repellent in h2
    probGoToh2[h1] = exp(-pow(distance[h1],2)/(2*pow(lambda,2))) * (1-static_cast<double>(combo2[h1])) ;
   
   //mosquitoes do not end up in the same household
   if (house2[h1] == house1[h1]) {probGoToh2[h1]=0.0;}
}

// sum of probabilities of moving to house h2 for each house h1
// zero out arrays
for (int h=0; h<30; h++) {sumProbGoToh2[h]=0;}
for (int h1=0; h1<(30*30); h1++) {
    m = (h1/30);
    sumProbGoToh2[m] = sumProbGoToh2[m] + probGoToh2[h1]; 
	}

// scale probabilities of moving to house h2 to one (conditional on coming from h1)
for (int h1=0; h1<(30*30); h1++) {
    m = (h1/30);
	probGoToh2[h1] = probGoToh2[h1]/sumProbGoToh2[m];
	}
//absolute probability of h1 to h2
for (int h1=0; h1<(30*30); h1++) {
    incomingh2[h1] = probGoToh2[h1] * (outgoingh1[h1]);
    }

// summarize to one record per house
// get outgoing per house (h) and incoming per house (h)
// zero out array
  for (int h=0; h<30; h++) {
      incomingh[h]=0.0;
      outgoingh[h]=0.0;
  }
  for (int h1=0; h1<(30*30); h1++) {
     m = (h1 % 30);  
     incomingh[m] = incomingh[m] + incomingh2[h1] ; 
     m = h1/30;
     outgoingh[m] = outgoingh[m] + (outgoingh1[h1]/30);
 }

// get predicted probability (missings set to 0 since data will be zero if missing)
for (int h=0; h<30; h++) {
 	//calculate the predicted probabilities for each house
    predph[h] = (pr[h] - outgoingh[h] + (incomingh[h] * hseFactor)) * notmissing[h];
    }

//to avoid missing values in the denominator 
for (int h=0; h<30; h++) {
    predph[h]=predph[h]+0.000001;
}

// scale predph to take account of missings === and also account for the mosquitoes ending up elsewhere
double sumpredph = 0.0;
for (int h=0; h<30; h++) {
    sumpredph = sumpredph + predph[h];
}

//calculate the adjusted predicted probability
for (int h=0; h<30; h++) { 
    adjpredph[h] = predph[h]/sumpredph;
       }

  // multinomial log likelihood
  for (int h=0; h<31; h++) {
    //ignore if adjpredprob is zero (missings)
    if (adjpredpha[h]>0) {
           loglik = loglik + (obsmosq[h + (d*30)] *log(adjpredph[h]));

	 }

}

} // close loop for d

// write outputs to a results file
ofstream outfile("results.txt");

 // baseline distributions of mosquitoes (pr) and simulated numbers of mosquitoes per house (obsmosq) excluding those diverted elsewhere
for (int h=0; h<30; h++) {
	outfile << h << " " << pr[h] << " " << obsmosq[h] <<" " << endl;
  }

 // simulated numbers of mosquitoes and the adjusted predicted probability accounting for mosquitoes diverted elsewhere
for (int h=0; h<31; h++) {
	outfile << h << " " << obsmosq[h] << "  "<< adjpredph[h] << "  "<< endl;
  }

 outfile << "repFactor " << spatialRepFactor << " " << endl; // proportion of mosquitoes repelled from households using repellents
 outfile << "lambda " << lambda << " " << endl; // mean distance moved between houses 
 outfile << "hseFactor " << hseFactor << " " << endl; // proportion of mosquitoes diverted elsewhere
 outfile << "loglik " << loglik << " " << endl; // the loglikelihood
 // outfile << "dissapear " << dissapear << " " << endl; // total number of mosquitoes diverted elsewhere
 // outfile << "ntotal " << ntotal << " " << endl; // total number of mosquitoes that exist

 outfile.close();

 return(-loglik);

}



