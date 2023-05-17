#include <iostream>
#include <cstdlib>
#include <vector>
#include <set>
#include <random>
#include <algorithm>
#include <fstream>
#include <chrono>



// ten main10_GR jest przystosowany do wczytywania listy polaczen .edgelist oraz listy wezlow wraz z polozeniami na kole. Odleglosci w przestrzeni ukrytej sa liczone przez GR na biezaco.
// main11_GR jest przystosowany do wczytywania listy polaczen i na jej podstawie stwierdza on liczbe N, sluzy on do wczytywania wylacznie sieci zlozonych z pojedynczego komponentu
// main16_GR przystosowuje do dzialania nie tylko na spojnym komponencie (chodzi o to ze wezly o stopniu 0 nie byly wspomniane na liscie polaczen), ale takze na pelnej sieci gdzie moga wystapic
// ... niewspomniane w ogole na liscie wezly. W tym celu podaje dodatkowo jako argument liczbe wezlow sieci N.
// main18_GR zmiemian dystans przechowywany w connection matrix z double na float
// main19_GR to samo co main18, z jedyna zmiana ze dodaje wypisanie ostrzezenia, co w tej chwili jest amalizowane: cala siec czy najwiekszy klaster
// main20_GR do pliku wynikowego gr summary dodaje kolumne z nazwa danej instancji sieci
// main21_GR w tym mainie zmieniam sposob wczytywania plikow tak, zeby liczone byly tylko pliki .edgelist, tj. zeby w folderze mogly znajdowac sie takze inne pliki
// main22_GR dodaje mozliwosc trzeciego argumentu - numer folderu

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



int greedy_routing(const unsigned N, const vector<vector<pair<bool,float>>>& connectivity_matrix, const int node_start, const int node_target, const double R)
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



unsigned findMax(set<unsigned> my_set)
{
 
    // Get the maximum element
    unsigned max_element;
    if (!my_set.empty())
        max_element = *(my_set.rbegin());
 
    // return the maximum element
    return max_element;
}






int main(int argc, char *argv[])
{

	if(argc != 1 and argc !=2 and argc!=3 and argc!=4)
	{
		cerr << "bad number of arguments" << endl;
		
		cerr << "Arguments shall be:" <<endl;
		cerr << "1. numer folderu _x, 2. liczba powtorzen GR dla każdej instancji sieci, 3. wielkosc sieci N " << endl; 
		return 1;
	}



	string numer_folderu = string("_1");		
	if(argc >= 2)
		numer_folderu = string(argv[1]);


	int max_gr = 10e4;;
	if(argc >= 3)
		max_gr = atoi(argv[2]);
	

	vector<double> wyniki_gr_kontener_instacje;
	vector<string> wyniki_gr_kontener_instacje_nazwy;


	namespace fs = std::filesystem;

	int fileCount = 0;
	
	unsigned N_component;
	//double alfa;
	//double gamma;

	std::string path = "./data_topologie" + numer_folderu; // uwaga, zmieniajac te sciezce pamietaj zmienic "extract erase", linie 373+
	
	
	//zliczmy pliki z tego folderu
	auto dirIter = fs::directory_iterator(path);

	/*
	int overallFileCount = std::count_if(
	    begin(dirIter),
	    end(dirIter),
	    [](auto& entry) { return entry.is_regular_file(); }
	);
	
	overallFileCount /= 2; // uwzgledniam ze plikow jest 2 razy wiecej
	*/
	
	int overallFileCount = 0;
	
	for (const auto & entry : fs::directory_iterator(path))
	{
		if( entry.path().extension() == ".edgelist" )
			overallFileCount++;
	}
	
	
	for (const auto & entry : fs::directory_iterator(path))
	{
		
		if( entry.path().extension() == ".txt" ) // .txt to pliki z wezlami, pozniej przeanalizujemy odpowiedni plik z wezlami
    			continue;
    			
    		if( entry.path().extension() != ".edgelist" ) // .edgelist to pliki z lista polaczen, w tej chwili tylko je chcemy 
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
		
		
		//input_file >> N;
		//input_file >> gamma;
		//input_file >> alfa;
		
		
		
		
		
		vector<unsigned> skad_, dokad_, wspomniane_wezly_;
		
		unsigned skad;
		unsigned dokad;
		
		double odleglosc_unused;
		
		while (input_file >> skad)
		{
			input_file >> dokad;
			input_file >> odleglosc_unused;
			
			skad_.push_back(skad);
			dokad_.push_back(dokad);
			
			wspomniane_wezly_.push_back(skad);
			wspomniane_wezly_.push_back(dokad);
		}
		
		
		input_file.close();
		
		
		//cout << "zamknalem plik" << endl;
		
		// liczba N_component jako liczba unikalnych elementow w zbiorze wspomnianych wezlow
		
		set<unsigned> zbior_numerow_wezlow( wspomniane_wezly_.begin(), wspomniane_wezly_.end() );
		
		N_component = zbior_numerow_wezlow.size();
		
		int max_num_wezla = findMax(zbior_numerow_wezlow);
		
		
		wspomniane_wezly_.clear();
		
		
		int N_podane;
		// jesli podano 3. argument tj rozmiar sieci - wowczas chcemy dzialac na calej sieci - wymusmy dzialanie na calej sieci
		if(argc >= 4)
		{
			N_podane = atoi(argv[3]);
			max_num_wezla = N_podane-1;
			
			wspomniane_wezly_ = vector<unsigned>(N_podane);
			
			for(int i=0; i<N_podane; i++)
				wspomniane_wezly_.at(i) = i;
				
			zbior_numerow_wezlow = set<unsigned>(wspomniane_wezly_.begin(), wspomniane_wezly_.end() );
			//N_component = zbior_numerow_wezlow.size();
			N_component = N_podane;
			
			cout << "Analizuje cala siec." << endl;
		}
		else
		{
			cout << "Analizuje tylko najwiekszy komponent." << endl;
		}
		
		
		vector<vector<pair<bool,float>>> connectivity_matrix(max_num_wezla+1);

		for(int i=0; i<connectivity_matrix.size(); i++)
			connectivity_matrix.at(i) = vector<pair<bool,float>>(max_num_wezla+1, make_pair(false, 0.));
		
		
		//cout << "zamknalem plik 2" << endl;
		
		//unsigned skad;
		//unsigned dokad;
		//bool polaczenie;
		//double odleglosc;
		
		
		/*
		// przebiegam po wektorach skad_ i dokad_ i 
		
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
			
			
			
			
			
		}
		/**/
		
		
		
		
		cout << "Wczytuje dane o polozeniach wezlow i obliczam odleglosci ukryte" << endl;
		
		
		vector<Node> nodes_;
		
		
		// extract current number of network
		
		string extract = entry.path();
		extract.erase(0, 27);
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
		
		//cout << "zamknalem plik 2.5" << endl;
		
		// zakladam ze lista nodes jest posortowana rosnaco !
		
		while (input_file_nodes >> numer_wezla)
		{
			input_file_nodes >> theta;
			
			
			Node node;
			node.theta = theta;
			
			nodes_.push_back( node );
				
		}
		
		input_file_nodes.close();
		
		//cout << "zamknalem plik 3" << endl;
		
		// teraz mam informacje o wektorze skad_, dokad_, mam nowe N_component (mniejsze), mam thete dla kazdego oryginalnego numeru wezla
		
		// przebiegam po wektorach skad_ i dokad_ i ustawiam polaczenia w macierzy polaczen - bede musial pamietac zeby losowac numery wezlow tylko z listy dostepnych wezlow
		
		for(int i=0; i<skad_.size(); i++)
		{
			unsigned skad = skad_.at(i);
			unsigned dokad = dokad_.at(i);
			
			//cout << skad << "  " << dokad << "  " << "  " << max_num_wezla <<   endl;
			//getchar();
			
			connectivity_matrix.at(skad).at(dokad).first = true;
			//connectivity_matrix.at(skad).at(dokad).second = odleglosc;
			
			connectivity_matrix.at(dokad).at(skad).first = true;
			//connectivity_matrix.at(dokad).at(skad).second = odleglosc;
			
			
		}
		
		
		//cout << "zamknalem plik 4" << endl;
		
		
		// teraz policzmy wszystkie odleglosci miedzy wezlowe -to musi byc zrobione juz po utworzeniu mneijszej macierzy 
		
		int N = nodes_.size(); // wielkosc oryginalnej sieci
		double R = N/(2.*PI);
		
		for(int i=0; i<nodes_.size(); i++)
			for(int j=i; j<nodes_.size(); j++) 
			{
				if(i==j) continue;
				
				// jesli numer wezla jest wiekszy niz numer wezla maksymalny jaki wystepuje w komponencie, to omijamy go, bo i tak nie ma gdzie go wpisac
				// (macierz ma rozmiar max_num_wezla x max_num_wezla)
				// czesc odleglosci (dla tych numerow wezlow ktore nie wystepuja w komponencie i tak sie nie przyda)
				if(i>max_num_wezla or j>max_num_wezla) continue;

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

		// teraz zamiast losowac dwoch numerow od 0 do N-1 musze losowac dwa numery sposrod dostepnych w tym komponencie numerow wezlow
		// numery wezlow dostepne sa w zbiorze zbior_numerow_wezlow
		// wylosuje wiec dwa randomowe miejsca w zbiorze, po czym wezme numery wezlow z tych miejsc


		for(int numer_gr = 0; numer_gr < max_gr; numer_gr++)
		{

			//cout << numer_gr << "  " << max_gr << endl;


			int wezel_a;
			int wezel_b;

			// wylosuj pare wezlow
			do
			{
				wezel_a = distribution(generator) % N_component ;
				wezel_b = distribution(generator) % N_component ;	
			}
			while(wezel_b == wezel_a);
			
			//przetlumacz miejsca w zbiorze na numery wezlow
			wezel_a = *next(zbior_numerow_wezlow.begin(), wezel_a);
			wezel_b = *next(zbior_numerow_wezlow.begin(), wezel_b);
			
			// uruchamiam greedy routing !

			int gr_result = greedy_routing(max_num_wezla+1, connectivity_matrix, wezel_a, wezel_b, R);
			
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
		
		wyniki_gr_kontener_instacje.push_back( sredni_success_ratio );
		wyniki_gr_kontener_instacje_nazwy.push_back( "network_" + extract  );

	} // end of for over files
	
	
	
	// obliczam srednie success ratio dla danego N, gamma i alfa
	
	double sredni_s_r = 0.0;

	for(auto& w:wyniki_gr_kontener_instacje)
		if(w > 0)
			sredni_s_r += w;

	sredni_s_r /= fileCount;

	//cout << endl << "sredni p_s dla N=" << N << ", gamma=" << gamma << ", i alfa=" << alfa << " jest " << sredni_s_r << endl;
	cout << endl << "sredni p_s przy tych N, alfa i gamma (nie wiem jakich konkretnie) jest " << sredni_s_r << endl;
	/**/
	
	

	//dla kazdego N, gamma, alfa : dodatkowy plik zbiorczy jak zachowywalo sie ps dla kolejnych realizacji sieci po uruchomieniu greedy routing
	
	//string filename_grsummary = "grsummary_N=unkn" + to_string_with_precision(N_component) + "_gamma=unkn" + "_alfa=unkn" + ".txt";
	string filename_grsummary = "grsummary_N=unkn_gamma=unkn_alfa=unkn.txt";
	ofstream output_grsummary(filename_grsummary.c_str());
		
	for(int i=0; i<fileCount; i++)
	{
		output_grsummary << wyniki_gr_kontener_instacje_nazwy.at(i) << "\t" << wyniki_gr_kontener_instacje.at(i) << endl; 
	
	}
		
	output_grsummary.close();
	/**/


	






	return 0;
}
