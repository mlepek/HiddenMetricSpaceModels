#include <iostream>
#include <cstdlib>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>


using namespace std;


struct Node
{
	double theta = 0.0;
	double kappa = 0.0;
	int k = 0;
};


double PI = 3.141592653589793238;


inline double ro_dist(double gamma, double avg_k, double kappa_0, double kappa)
{
	return (gamma - 1.) * pow(kappa_0, gamma-1.) * pow(kappa, -gamma);
}


double d_theta_theta(double theta_1, double theta_2, double R)
{
	// lets make sure that theta_2 is greater than theta_1
	/*if(theta_1 > theta_2)
		swap(theta_1, theta_2);

	double distance = (theta_2 - theta_1) * R;

	double alternative_dist = 2*PI-theta_2 + theta_1;

	if(alternative_dist < distance)
		distance = alternative_dist;
	*/

	// wedle ksiazki
	double distance = R * (PI - abs(PI - abs(theta_1 - theta_2)));

	return distance;
}








int main(int argc, char *argv[])
{

	if(argc != 4)
	{
		cerr << "bad number of arguments" << endl;
		return 1;
	}

	double gamma = atof(argv[1]);
	double alfa = atof(argv[2]);
	unsigned N = atoi(argv[3]);

	double avg_k = 6.0;

	// creating vector of nodes
	vector<Node> nodes_(N);

	// initialization of the random number generator
	default_random_engine generator;
	int maxInt = 2147483647;
	uniform_int_distribution<int> distribution(0, maxInt-1);

	// for random seed
	generator.seed( time(NULL) );

	// filling vector of nodes with random theta

	for(auto& n : nodes_)
	{
		// getting random double from 0 - 2 pi
		n.theta = double( distribution(generator) ) / maxInt * 2. * PI;
	}

	// filling vector of nodes with random kappa from ro distribution

	for(auto& n : nodes_)
	{
		//

		double kappa;

		while(true)
		{
			double kappa_0 = (gamma - 2.) * avg_k / (gamma - 1.);
			
			kappa = double( distribution(generator) ) / maxInt * 1.e4 + kappa_0; // random double from kappa_0 to 10 000 + kappa_0

			//cout << kappa << "  ";

			double max_y = ro_dist(gamma, avg_k, kappa_0, kappa_0); //// ????????

			//cout << max_y << "  ";

			double y = double( distribution(generator) ) / maxInt * max_y;
			//double y = double( distribution(generator) ) / maxInt * 100;



			double ro = ro_dist(gamma, avg_k, kappa_0, kappa);

			//cout << y << "  " << ro << endl;

			//cout << kappa_0 << "  " << gamma-1. << "  " <<  pow(kappa_0, gamma-1.) << "  " <<  pow(kappa, -gamma) << "  " << ro << "  " << y << endl;
			//cout << y << "  " << ro << endl;

			if( y <= ro )
				break;
		}

		n.kappa = kappa;
	}
	

	double srednie_kappa = 0.0;

	for(auto& n : nodes_)
		srednie_kappa += n.kappa;

	srednie_kappa /= N;

	cout << "srednie kappa jest " << srednie_kappa << ", a powinno byc okolo " << avg_k << endl;


	// -----------------------------------------------------


	vector<vector<bool>> connectivity_matrix(N);

	for(int i=0; i<connectivity_matrix.size(); i++)
		connectivity_matrix.at(i) = vector<bool>(N, false); // we will use only the upper-right half


	//double mju = (alfa-1.)/avg_k;
	double mju = alfa / (2*PI*avg_k) * sin(PI/alfa); // wg ksiazki

	double R = N/(2.*PI);


	for(int i=0; i<nodes_.size(); i++)
		//for(int j=i+1; j<nodes_.size(); j++) // we treat any pair only once and do not treat i=j
		for(int j=i; j<nodes_.size(); j++) 
		{
			if(i==j) continue;

			double dist = d_theta_theta(nodes_.at(i).theta, nodes_.at(j).theta, R);


			// probability of connection

			//double argument = 1. + dist / (mju * nodes_.at(i).kappa * nodes_.at(j).kappa);
			//double r = pow(argument, -alfa);

			// wg ksiazki
			double r = 1. / ( 1. + pow( dist / (mju * nodes_.at(i).kappa * nodes_.at(j).kappa), alfa ) );

			double random_float = double( distribution(generator) ) / maxInt;

			if( random_float < r )
			{
				// para i oraz j ma zostac polaczona
				connectivity_matrix.at(i).at(j) = true;
				connectivity_matrix.at(j).at(i) = true;
			}

		}	 




	// sprawdzic czy rzeczywiscie sredni stopien jest <k> = 6

	for(int i=0; i<N; i++)
	{

		int k = 0;

		for(int j=0; j<N; j++)
		{
			// sumowac polaczenia danego wierzcholka
			if ( connectivity_matrix.at(i).at(j) == true )
				k++;
		}

		nodes_.at(i).k = k;
	}

	double srednie_k = 0.0;

	for(auto& n : nodes_)
		srednie_k += n.k;

	srednie_k /= double(N);

	cout << "srednie k rzeczywiste jest " << srednie_k << endl;


	// histogram

	vector<int> histogram(N);

	for(auto& n : nodes_)
	{
		histogram.at( n.k )++;
	}

	//for(auto& h : histogram)
	//	cout << h << "  " ;
	cout << endl;



	// sprawdzic czy rezczywiscie rozklad k jest jak k^-gamma


	ofstream output_file("hist_out.txt");
	for(auto& h : histogram)
		output_file << h << endl;
	output_file.close(); 


	// warto by policzyc gamma tak jak jest u MEJ w "Zipf Law"

	// losowanie z "unbounded distribution" lepsze niz za pomocą punktu z 2-wymiarowej powierzchni

	// zastosowanie rownan z cutoff i przeciwko finite size effects

	// zasymulowanie greedy routing dla duzych gamma, 2.8 - 3.0 i duzych alfa 3 - 5



	//for(auto n : nodes_)
	//	cout << n.kappa << endl;

	/*
	double srednia_kappa = 0.0;

	for(auto& n:nodes_)
	{
		srednia_kappa += n.kappa;
	}

	srednia_kappa /= nodes_.size();

	cout << srednia_kappa << endl;
	*/







	return 0;
}