#include <iostream>
#include <cstdlib>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
#include <chrono>



// ten main10_GR jest przystosowany do wczytywania listy polaczen .edgelist oraz listy wezlow wraz z polozeniami na kole. Odleglosci w przestrzeni ukrytej sa liczone przez GR na biezaco.


#include <sstream>
#include <filesystem>


template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 1)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}


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


// wersja z poprawkami przeciwko finite-size effects wg ksiazki (against finite size effects, "afse")

inline double ro_dist_afse(double gamma, double avg_k, double kappa_0, double kappa, double kappa_c)
{
	return (gamma - 1.) * pow(kappa_0, gamma-1.) * pow(kappa, -gamma) / (1. - pow(kappa_c/kappa_0, 1.-gamma) );
}


inline double d_theta_theta(double theta_1, double theta_2, double R)
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
	return R * (PI - abs(PI - abs(theta_1 - theta_2)));
}



int greedy_routing(const unsigned N, const vector<vector<pair<bool,double>>>& connectivity_matrix, const int node_start, const int node_target, const double R)
{

	int current_node = node_start;

	vector<bool> nodes_visited_(N, false);

	nodes_visited_.at(node_start) = true;

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

	if(argc != 1 and argc !=2)
	{
		cerr << "bad number of arguments" << endl;
		
		cerr << "Arguments shall be:" <<endl;
		cerr << "liczba powtorzen GR dla każdej instancji sieci" << endl; 
		return 1;
	}







	int max_gr;
		if(argc == 2)
			max_gr = atoi(argv[1]);
		else
			max_gr = 10e4;

	vector<double> wyniki_gr_kontener_instacje;


	namespace fs = std::filesystem;

	int fileCount = 0;
	
	unsigned N;
	double alfa;
	double gamma;

	std::string path = "./data_topologie"; // uwaga, zmieniajac te sciezce pamietaj zmienic "extract erase", linie 271-275
	
	
	//zliczmy pliki z tego folderu
	auto dirIter = fs::directory_iterator(path);

	int overallFileCount = std::count_if(
	    begin(dirIter),
	    end(dirIter),
	    [](auto& entry) { return entry.is_regular_file(); }
	);
	
	overallFileCount /= 2; // uwzgledniam ze plikow jest 2 razy wiecej
	
	
	for (const auto & entry : fs::directory_iterator(path))
	{
		
		if( entry.path().extension() == ".txt" ) // patrzymy na plik z wezlami
    			continue;
		
	
		fileCount++;
		//std::cout << entry.path() << std::endl;
		
		cout << endl << "plik " << fileCount << "/" << overallFileCount << endl;
		
		std::cout << "analizuje " <<  entry.path() << std::endl;
		
		ifstream input_file( entry.path().c_str() );
		
		
		if (!input_file.is_open())
		{
			cerr << "Could not open the file - '" << entry.path() << "'" << endl;
			return EXIT_FAILURE;
		}
		
		
		input_file >> N;
		input_file >> gamma;
		input_file >> alfa;
		
		
		vector<vector<pair<bool,double>>> connectivity_matrix(N);

		for(int i=0; i<connectivity_matrix.size(); i++)
			connectivity_matrix.at(i) = vector<pair<bool,double>>(N, make_pair(false, 0.));
		
		
		unsigned skad;
		unsigned dokad;
		//bool polaczenie;
		double odleglosc;
		
		while (input_file >> skad)
		{
			input_file >> dokad;
			//input_file >> polaczenie;
			input_file >> odleglosc;
			
			//cout << skad << "  " << dokad << "  " << odleglosc << endl;
			//getchar();
			
			connectivity_matrix.at(skad).at(dokad).first = true;
			connectivity_matrix.at(skad).at(dokad).second = odleglosc;
			
			connectivity_matrix.at(dokad).at(skad).first = true;
			connectivity_matrix.at(dokad).at(skad).second = odleglosc;
			
			
			
			/**/
			
		}
		
		
		input_file.close();
		
		
		cout << "Wczytuje dane o polozeniach wezlow i obliczam odleglosci ukryte" << endl;
		
		
		vector<Node> nodes_(N);
		
		
		// extract current number of network
		
		string extract = entry.path();
		extract.erase(0, 25);
		string to_find = "_";
		size_t found = extract.find(to_find);
		extract = extract.substr(0, found);
		
		//cout << extract << endl;
		//getchar();
		
		// teraz wczytajmy plik z wezlami i ich poloizeniami na kole
		
		string filename_nodes = path + "/network_" + extract + "_nodes.txt";
		ifstream input_file_nodes( filename_nodes.c_str() );
		
		
		if (!input_file_nodes.is_open())
		{
			cerr << "Could not open the file - '" << filename_nodes << "'" << endl;
			return EXIT_FAILURE;
		}
		
		unsigned numer_wezla;
		double theta;
		
		while (input_file_nodes >> numer_wezla)
		{
			input_file_nodes >> theta;
			
			nodes_.at(numer_wezla).theta = theta;	
		}
		
		input_file_nodes.close();
		
		
		
		// teraz policzmy wszystkie odleglosci miedzy wezlowe
		
		double R = N/(2.*PI);
		
		for(int i=0; i<nodes_.size(); i++)
			for(int j=i; j<nodes_.size(); j++) 
			{
				if(i==j) continue;

				double dist = d_theta_theta(nodes_.at(i).theta, nodes_.at(j).theta, R);

				connectivity_matrix.at(j).at(i).second = dist;
				connectivity_matrix.at(i).at(j).second = dist;
			}	 

		// w ten sposob mamy przygotowana connectivity matrix do greedy routing: polaczenia sa tylko tam gdzie trzeba, a adleglosci policzone sa dla wszystkich wezlow

		
		
		
		
		
		/*
		for(int i=0; i<N; i++)
			{
			for(int j=0; j<N; j++)
			{
				 
					cout << connectivity_matrix.at(i).at(j).second << " " ; 
					
			
			}
			cout << endl;
			}
		getchar();
		/**/
		
		/*
		for(int i=0; i<N; i++)
		{
			for(int j=i+1; j<N; j++)
			{
				if( connectivity_matrix.at(i).at(j).first == true )
					cout << i << " " << j << " " << connectivity_matrix.at(i).at(j).second << endl; 
			}
		}
		
		getchar();
		/**/
		
		
		cout << "uruchamiam greedy routing" << endl;
		
		// initialization of the random number generator
		default_random_engine generator;
		int maxInt = 2147483647;
		uniform_int_distribution<int> distribution(0, maxInt-1);

		// for random seed
		auto seed = static_cast<long unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
		generator.seed( seed );
		
		
		//double R = N/(2.*PI);
		
		
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
			

			int gr_result = greedy_routing(N, connectivity_matrix, wezel_a, wezel_b, R);
			
			//cout << N << endl;

			//cout << "instancja nr " << fileCount << ", powtorzenie gr " << numer_gr << ", hop length: " << gr_result << endl;

			wyniki_gr.push_back(gr_result);

		}
		
		
		// obliczam sredni success ratio dla danej instacji sieci

		double sredni_success_ratio = 0.0;

		for(auto& w:wyniki_gr)
			if(w > 0)
				sredni_success_ratio++;

		//cout << "ssr " << sredni_success_ratio << "  , max_gr " << max_gr << endl;

		sredni_success_ratio /= (max_gr-1);

		//cout << "ssr " << sredni_success_ratio << "  , max_gr " << max_gr << endl;

		//getchar();
		cout << "plik " << entry.path() << ", sredni s.r.: " << sredni_success_ratio << endl;
		
		wyniki_gr_kontener_instacje.push_back(sredni_success_ratio);

	} // end of for over files
	
	
	
	// obliczam srednie success ratio dla danego N, gamma i alfa
	
	double sredni_s_r = 0.0;

	for(auto& w:wyniki_gr_kontener_instacje)
		if(w > 0)
			sredni_s_r += w;

	sredni_s_r /= fileCount;

	cout << endl << "sredni p_s dla N=" << N << ", gamma=" << gamma << ", i alfa=" << alfa << " jest " << sredni_s_r << endl;
	/**/
	
	

	//dla kazdego N, gamma, alfa : dodatkowy plik zbiorczy jak zachowywalo sie ps dla kolejnych realizacji sieci po uruchomieniu greedy routing
	
	string filename_grsummary = "grsummary_N=" + to_string_with_precision(N) + "_gamma=" + to_string_with_precision(gamma) + "_alfa=" + to_string_with_precision(alfa) + ".txt";
	ofstream output_grsummary(filename_grsummary.c_str());
		
	for(int i=0; i<fileCount; i++)
	{
		output_grsummary << i << " " << wyniki_gr_kontener_instacje.at(i) << endl; 
	
	}
		
	output_grsummary.close();
	/**/


	






	return 0;
}
