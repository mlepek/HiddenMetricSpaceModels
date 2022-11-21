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



int greedy_routing(const vector<Node>& nodes_, const vector<vector<pair<bool,double>>>& connectivity_matrix, const int node_start, const int node_target, const double R)
{

	int current_node = node_start;

	vector<bool> nodes_visited_(nodes_.size(), false);

	nodes_visited_.at(node_start) = true;

	//cout << nodes_visited_.size() << " " << nodes_.size() << "   " << node_start << " " << node_target << endl;

	int hop_length = 0;

	while(true)
	{
		if( current_node == node_target )
			break;

		// stwierdzmy, ktory wierzcholek z dostepnych wierzcholkow ma najkrotszy dystans ukryty do wierzcholka target

		int index_of_neighbour_closest_to_target = -1;
		//cout << index_of_neighbour_closest_to_target << endl;

		double minimal_distance = R*100000000; // na pewno wiekszy od jakiekolwiek mozliwego dystansu

		// przejazd po wierzcholkach sasiadach current node'a
		for(int i=0; i<connectivity_matrix.at(current_node).size(); i++ )
		{
			if(connectivity_matrix.at(current_node).at(i).first) // jesli w ogole jest to sasiad 
			{
				// sprawdz jaka ma odleglosc ukryta do target node

				if(connectivity_matrix.at(node_target).at(i).second < minimal_distance)
				{

					// zapamietaj ow dystans i numer wierzcholka
					minimal_distance = connectivity_matrix.at(node_target).at(i).second;
					index_of_neighbour_closest_to_target = i;
				}
			}
		}

		if(index_of_neighbour_closest_to_target==-1) // tzn ze w ogole nie znaleziono sasiada blizszego targetowi - nie wiem czy w ogole powinno dochodzic do takiej sytuacji
			index_of_neighbour_closest_to_target = current_node;



		// w tej chwili mam wybrany wierzcholek najblizszy w przestrzeni ukrytej targetowi

		current_node = index_of_neighbour_closest_to_target;

		// sprawdzmy czy juz bylismy w tym wierzcholku
		if(nodes_visited_.at(current_node) == true)
			return -1; // greedy routing wpadl w petle
		else
		{
			nodes_visited_.at(current_node) = true;
			hop_length++;
		}

	}


	return hop_length;
}









int main(int argc, char *argv[])
{

	if(argc != 5)
	{
		cerr << "bad number of arguments" << endl;
		return 1;
	}

	double gamma = atof(argv[1]);
	double alfa = atof(argv[2]);
	unsigned N = atoi(argv[3]);




	int max_iteracji = atoi(argv[4]); // liczba instancji sieci do wygenerowania


	//double gamma = 3.0;
	//double alfa = 5.0;
	//unsigned N = 1000;

	double avg_k = 6.0;

	vector<double> wyniki_gr_kontener_instacje;


	for(int iteracje=0; iteracje<max_iteracji; iteracje++)
	{


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

		cout << endl << "instancja " << iteracje <<  ", srednie kappa jest " << srednie_kappa << ", a powinno byc okolo " << avg_k << endl;


		// -----------------------------------------------------


		vector<vector<pair<bool,double>>> connectivity_matrix(N);

		for(int i=0; i<connectivity_matrix.size(); i++)
			connectivity_matrix.at(i) = vector<pair<bool,double>>(N, make_pair(false, 0.)); // we shall maybe use only the upper-right half


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
					connectivity_matrix.at(i).at(j).first = true;
					connectivity_matrix.at(j).at(i).first = true;
				}

				// odleglosc w przestrzeni ukrytej zapisywana jest zawsze
				connectivity_matrix.at(j).at(i).second = dist;
				connectivity_matrix.at(i).at(j).second = dist;

			}	 




		// sprawdzic czy rzeczywiscie sredni stopien jest <k> = 6

		for(int i=0; i<N; i++)
		{

			int k = 0;

			for(int j=0; j<N; j++)
			{
				// sumowac polaczenia danego wierzcholka
				if ( connectivity_matrix.at(i).at(j).first == true )
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
		//cout << endl;



		// sprawdzic czy rezczywiscie rozklad k jest jak k^-gamma

		/*
		ofstream output_file("hist_out.txt");
		for(auto& h : histogram)
			output_file << h << endl;
		output_file.close(); 
		*/


		int max_gr = 10e4; // w pracy 10e6


		vector<int> wyniki_gr;


		for(int numer_gr = 0; numer_gr < max_gr; numer_gr++)
		{

			//cout << numer_gr << "  " << max_gr << endl;


			int wezel_a;
			int wezel_b;

			// wylosuj pare wezlow
			do
			{
				wezel_a = distribution(generator) % N ;
				wezel_b = distribution(generator) % N ;	
			}
			while(wezel_b == wezel_a);
			

			int gr_result = greedy_routing(nodes_, connectivity_matrix, wezel_a, wezel_b, R);

			//cout << "instancja nr " << iteracje << ", powtorzenie gr " << numer_gr << ", hop length: " << gr_result << endl;

			wyniki_gr.push_back(gr_result);

		}

		//cout << "instancja nr " << iteracje << ", rozmiar wyniki_gr_ : " << wyniki_gr.size() << endl;

		//for(auto& w:wyniki_gr)
		//	cout << w << " ";
		//cout << endl;
		//getchar();

		// obliczam sredni success ratio dla danej instacji sieci

		double sredni_success_ratio = 0.0;

		for(auto& w:wyniki_gr)
			if(w > 0)
				sredni_success_ratio++;

		//cout << "ssr " << sredni_success_ratio << "  , max_gr " << max_gr << endl;

		sredni_success_ratio /= max_gr;

		//cout << "ssr " << sredni_success_ratio << "  , max_gr " << max_gr << endl;

		//getchar();
		cout << "instancja nr " << iteracje << ", sredni s.r.: " << sredni_success_ratio << endl;


		wyniki_gr_kontener_instacje.push_back(sredni_success_ratio);



	} ///// koniec for po instacjach sieci


	// obliczam srednie success ratio dla danego N, gamma i alfa

	double sredni_s_r = 0.0;

	for(auto& w:wyniki_gr_kontener_instacje)
		if(w > 0)
			sredni_s_r += w;

	sredni_s_r /= max_iteracji;

	cout << "sredni p_s dla N=" << N << ", gamma=" << gamma << ", i alfa=" << alfa << " jest " << sredni_s_r << endl;
	

	

	



	// zasymulowanie greedy routing dla duzych gamma, 2.8 - 3.0 i duzych alfa 3 - 5



	// warto by policzyc gamma tak jak jest u MEJ w "Zipf Law"

	// losowanie z "unbounded distribution" lepsze niz za pomocą punktu z 2-wymiarowej powierzchni

	// zastosowanie rownan z cutoff i przeciwko finite size effects

	



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