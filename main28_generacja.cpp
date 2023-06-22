#include <iostream>
#include <cstdlib>
#include <vector>
#include <list>
#include <random>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <cmath>

#include "cantor2.h"

// "Log" zmian:
// main7 pochodzi main5, w main7 wprowadzam rownania przeciwko "finite size effects" dla gamma bliskiego 2 od prawej strony
// main8 to podzial na dwa programy, generujacy sieci i wczytujace sieci z uruchomieniem greedy routing
// main9 tym razem testuje efektywna gamme (slope) histogramu stopni wezlow w wygenerowanej sieci
// main10 tym razem, podobnie jak w main8, przygotowuje dwa programy: generacja, generujacy siec i outputujacy liste linkow oraz liste wezlow wraz z polozeniami na kole
// ... oraz drgui program, GR, wczytujacy liste linkow oraz liste wezlow wraz z polozeniami na kole (GR wymaga wiedzy o wszystkich odleglosciach miedzywezlowych, takze niepolaczonych)
// main11 nie zapisuje N, alfa, gamma do pliku z lista polaczen
// main12 dodaje generacje zbioru cantora i mozliwosc generacji polozen katowych wezlow na odcinkach kola okreslonych zbiorem cantora
// main13  tej edycji dodaje do kodu sledzenie rozmiarow klastrow i przynaleznosci wezlow do danych klastrow w trakcie agregacji tak, jak niegdys w przypadku symulacji agregacji
// ...i symulacji regul odwrotnych, co ma na celu latwy wybor najwiekszego klastra z wygenerowanej sieci. W topology edgelist wypisuje tylko wezly z najw. klastra.
// main15 w tym main wylaczam zbior cantora i wypisuje stopnie wierzcholkow zarowno calej sieci jak i najwiekszego klastra, zeby pozniej w pythonie wyswietlac rozklady i obliczac wykladnik gamma w sensie MEJ Newmana (...py)
// main16 tutaj zmieniam zapis z zapisu wylacznie najwiekszego klastra na zapis calej sieci, co ma posluzyc zwalidowaniu programu main11_GR
// main17 rozwiazanie problemu z niemaleniem p_s dla malych gamm -> zmiana kappa_0, poprawka jest stosowana po obliczeni kappa_c a nie przed; zapisywana jest cala siec,
// ... w connection matrix przechowuje takze odleglosci
// main18 ponownie zaprzestaje przechowywania odleglosci w connection matrix, przelaczam na zapis tylko najwiekszych klastrow
// main19 to samo co main18, ale zapis calej sieci - dla Macka
// main22
// main23 w tym mainie zamieniam sposob przechowywania sieci z macierzy polaczen na liste polaczen, czyli liste wektorow sasiadow; dodaje argumnety na nazwe output folderu i na wypisywanie calej sieci lub tylko najwiekszego komponentu; zwalidowalem pod wzgledem wynikow pozytywnie, ale wydajnosc jest niska z powodu poczwornej petli przy add_connection
// main24 w tym mainie usprawniam petle kreacji polaczen, gdyz w main23 efektywnie byla tam petla 4-krotna za sprawa funkcji add_connection; tutaj add_connection zostala "wpieta" w kod glowny, jednak zysk czasowy wyszedl tylko okolo 20%, ale wyszedl. Trzeba jeszcze zrobic walidacje
// main25 w tym mainie wypisuje do pliku tylko wezly bez odleglosci - dla liczenia wymiarow fraktalnych za pomoca box_covering.ipynb
// main26 rozszerzam cantor.h i main do obslugi zbiorow cantora o dowolnym wymiarze fraktalnym, dodaje odpowiedni argument wejsciowy, dodatkowo argument wjesciowy na poziom zbioru Cantora
// main27 ...
// main28 na podstawie main26, dodaje przeliczanie polozen wezlow z przestrzeni S1 na przestrzen H2 wg rownan z pracy Garcia-Perez 2019 "Mercator..."


#include <sstream>

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 1)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}


using namespace std;

typedef std::vector<int> IVector;


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
	//double distance = R * (PI - abs(PI - abs(theta_1 - theta_2)));

	//return distance;
	
	return R * (PI - abs(PI - abs(theta_1 - theta_2)));
}

inline double delta_theta(double theta_1, double theta_2)
{	
	return PI - abs(PI - abs(theta_1 - theta_2));
}



//int greedy_routing(const vector<Node>& nodes_, const vector<vector<pair<bool,double>>>& connectivity_matrix, const int node_start, const int node_target, const double R)
int greedy_routing(const unsigned N, const vector<vector<pair<bool,float>>>& connectivity_matrix, const int node_start, const int node_target, const double R)
{

	int current_node = node_start;

	vector<bool> nodes_visited_(N, false);

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



void add_connection( list<vector<unsigned>>& lst_, unsigned where_i, unsigned where_j )
{	
	unsigned counter = 0;
	
	for(list<vector<unsigned>>::iterator it = lst_.begin(); it != lst_.end(); ++it)
	{
		if(counter == where_i)
		{
			// dodajmy w wektorze wezla i sasiada j
			(*it).push_back(where_j);
			
			//cout << (*it).back() << endl;
		
			break; // i wyjdzmy z fora	
		}
		
		counter++;
	}

	return;
}


bool check_connection(const list<vector<unsigned>>& lst_, unsigned where_i, unsigned where_j )
{	
	unsigned counter = 0;
	
	bool isconnection = false;
	
	for(list<vector<unsigned>>::const_iterator it = lst_.begin(); it != lst_.end(); ++it)
	{
		if(counter == where_i)
		{
			//for(auto a:(*it))
			//	cout << a;
			//cout << endl;
		
			// sprawdzmy czy j jest sasiadem i, tj. czy j zawiera sie w wektorze sasiadow i
			if( std::find( (*it).begin(), (*it).end(), where_j ) != (*it).end() )
			{	
				// jest sasiadem
				isconnection = true;
				//cout << "jest sasiadem" << endl;
				break;
			}
		}
		
		counter++;
	}
	
	return isconnection;
}



int get_degree(const list<vector<unsigned>>& lst_, unsigned node_i)
{
	unsigned counter = 0;
	int degree = 0;
	
	for(list<vector<unsigned>>::const_iterator it = lst_.begin(); it != lst_.end(); ++it)
	{
		if(counter == node_i)
		{
			// liczba elementow wektora dla wezla i jest stopniem wezla i - po prostu
			degree = (*it).size();
			break;
		}
		
		counter++;
	}
	

	return degree;
}






int main(int argc, char *argv[])
{

	if(argc != 7 and argc != 8 and argc != 9 and argc != 10)
	{
		cerr << "bad number of arguments" << endl;
		
		cerr << "Arguments shall be:" <<endl;
		cerr << "gamma \n alfa \n N \n liczba instancji sieci do wygenerowania \n srednie oczekiwane k \n fn / lc i.e. what to save, full network or largest comp. in edgelist \n folder zapisu \n dzielnik Cantora \n poziom zglebienia Cantora" << endl; 
		return 1;
	}

	double gamma = atof(argv[1]);
	double alfa = atof(argv[2]);
	unsigned N = atoi(argv[3]);




	int max_iteracji = atoi(argv[4]); // liczba instancji sieci do wygenerowania


	//double gamma = 3.0;
	//double alfa = 5.0;
	//unsigned N = 1000;

	//double avg_k = 6;
	double avg_k = atof(argv[5]);
	
	string what_to_save = string("fn");
	if(argc >= 7)
		what_to_save = string(argv[6]);
	
	string where_to_save = string("data_topologie_1");	
	if(argc >= 8)
		where_to_save = string(argv[7]);
		
	double dzielnik_Cantora = 3.;
	if(argc >= 9)
		dzielnik_Cantora = atof(argv[8]);
		
	int poziom_Cantora = 4.;
	if(argc == 10)
		poziom_Cantora = atoi(argv[9]);
		

	vector<double> wyniki_gr_kontener_instacje;


	for(int iteracje=0; iteracje<max_iteracji; iteracje++)
	{

		cout << endl << "instancja " << iteracje << endl;

		// creating vector of nodes
		vector<Node> nodes_(N);

		// initialization of the random number generator
		default_random_engine generator;
		int maxInt = 2147483647;
		uniform_int_distribution<int> distribution(0, maxInt-1);

		// for random seed
		auto seed = static_cast<long unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
		generator.seed( seed );

		// filling vector of nodes with random theta
		
		
		if( dzielnik_Cantora <= 2. )
		{
			// bez uzycie zbioru cantora
			
			for(auto& n : nodes_)
			{				
				// getting random double from 0 - 2 pi
				n.theta = double( distribution(generator) ) / maxInt * 2. * PI;
			}
		}
		else
		{
			// z uzyciem zbioru cantora
		
			// generate cantor set on a circle
			// od 0 do 1 znaczy odcinek od 0 do 2pi, n-ty stopien zaglebienia, dzielnik zbioru Cantora
			vector<double> cantor_set = buildCantorSetOnCircle(0, 1, poziom_Cantora, dzielnik_Cantora); 
			
			for(auto& n : nodes_)
			{		
				double random_theta;
				bool theta_found = false;
					
				// using Cantor set		
				do
				{
					// getting random double from 0 - 2 pi
					random_theta = double( distribution(generator) ) / maxInt * 2. * PI;
						
					////check if this number falls into the cantor set or not
					for(int i=0; i<cantor_set.size(); i+=2)
					{
						if( random_theta >= cantor_set.at(i) and random_theta <= cantor_set.at(i+1) )
						{
							theta_found = true;
							break;
						}
					}
				}
				while( !theta_found );
				/**/
				
				n.theta = random_theta;
			}
		}
		
		
	
		
		//ofstream out("out.txt");
		//for(auto n : nodes_)
		//	out << n.theta << endl;
		//out.close();
			
		//getchar();

		// filling vector of nodes with random kappa from ro distribution
		
		double kappa_0 = (gamma - 2.) * avg_k / (gamma - 1.);
		//double kappa_0 =  (1. - 1./N) / (1. - pow(N, (2.-gamma)/(gamma-1.))) * (gamma - 2.) * avg_k / (gamma - 1.);
				
		double kappa_c = kappa_0 * pow(N, 1. / (gamma - 1.));
				

		for(auto& n : nodes_)
		{
			//

			double kappa;

			while(true)
			{
				
				
				//kappa = double( distribution(generator) ) / maxInt * 1.e4 + kappa_0; // random double from kappa_0 to 10 000 + kappa_0
				kappa = double( distribution(generator) ) / maxInt * (kappa_c - kappa_0) + kappa_0; // random double from kappa_0 to kappa_c

				//cout << kappa << "  ";

				double max_y = ro_dist(gamma, avg_k, kappa_0, kappa_0); //// mam nadzieje ze to nic nie psuje, ograniczam w ten sposob zakres osi pionowej losowania, wychodzac z zalozenia, ze szczyt rozkladu jest na samym poczatku

				//cout << max_y << "  ";

				double y = double( distribution(generator) ) / maxInt * max_y;
				//double y = double( distribution(generator) ) / maxInt * 100;



				// modyfikujemy teraz kappa_0 zgodnie ze wz. (2.8) z ksiazki, jakkolwiek nie jest na 100% jasne czy nie powinnismy zrobic tego juz na poczatku tj. przed obliczeniem kappa_c - ale chyba jednak raczej tutaj, bo jak bylo na poczatku to cos bylo zle
				kappa_0 =  (1. - 1./N) / (1. - pow(N, (2.-gamma)/(gamma-1.))) * (gamma - 2.) * avg_k / (gamma - 1.);
				
				

				//double ro = ro_dist(gamma, avg_k, kappa_0, kappa);
				// uwzglednienie rownan przeciwko finite size effects
				double ro = ro_dist_afse(gamma, avg_k, kappa_0, kappa, kappa_c);
				

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
		
		/*
		cout << "kreuje connectivity matrix" << endl;

		//vector<vector<pair<bool,float>>> connectivity_matrix(N);
		vector<vector<bool>> connectivity_matrix(N);

		for(int i=0; i<connectivity_matrix.size(); i++)
			//connectivity_matrix.at(i) = vector<pair<bool,float>>(N, make_pair(false, 0.)); // we shall maybe use only the upper-right half
			connectivity_matrix.at(i) = vector<bool>(N, false); // we shall maybe use only the upper-right half


		cout << "kreuje connectivity matrix - koniec " << endl;
		*/
		
		// zakladam liste polaczen
		list<vector<unsigned>> connectivity_list( N, vector<unsigned>(0) );
		
		
		
		//double mju = (alfa-1.)/avg_k;
		double mju = alfa / (2*PI*avg_k) * sin(PI/alfa); // wg ksiazki

		double R = N/(2.*PI);
		
		
		// ---------------- inicjalizacja elementow do sledzenia wielkosci klastrow ----------------
		
		
		// This vector stores labels (ordinal numbers) of mers (clusters).
		// Size of this vector decreases in time as the number of mers decreases.
		IVector labels_of_mers;

		// This vector contains masses (sizes) of mers. Its length decreases in time.
		IVector masses_of_mers;

		// This vector contains the information on: which cluster (mer) does this monomer belongs to.
		// This vector is of constant length.
		IVector membership_of_monomers;

		// "monodisperse initial conditions" as we would say in aggregation studies
		   
		// labels from 1 to N
		labels_of_mers = IVector(N);
		for (int i = 0; i<N; ++i)
			labels_of_mers[i] = i + 1;

		// all masses of size 1
		masses_of_mers = IVector(N, 1);

		// all are distinct clusters
		membership_of_monomers = IVector(N);
		for (int i = 0; i<N; ++i)
			membership_of_monomers[i] = i + 1;
		
		// -----------------------------------------------------------------------------------------


		list<vector<unsigned>>::iterator it_i = connectivity_list.begin();
		
		

		for(int i=0; i<nodes_.size(); i++, it_i++)
		{
			if ( i%1000 == 0 ) cout << " " << double(i)/N*100 << " % " << endl;
			
			
			list<vector<unsigned>>::iterator it_j = connectivity_list.begin();
			advance(it_j, i);
			
			//for(int j=i+1; j<nodes_.size(); j++) // we treat any pair only once and do not treat i=j
			for(int j=i; j<nodes_.size(); j++, it_j++) 
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
					/*
					// para i oraz j ma zostac polaczona
					//connectivity_matrix.at(i).at(j).first = true;
					connectivity_matrix.at(i).at(j) = true;
					//connectivity_matrix.at(j).at(i).first = true;
					connectivity_matrix.at(j).at(i) = true;
					*/
					//add_connection( connectivity_list, i, j );
					//add_connection( connectivity_list, j, i );
					
					// dodajmy w wektorze wezla i sasiada j
					(*it_i).push_back(j);
					
					// oraz na odwrot
					(*it_j).push_back(i);
					
					
					
					// w tym miejscu musi nastapic laczenie wezlow w sensie agregacji tj., zarzadzenie obecnymi wielkosciami klastrow, istnieniem klastrow
					// ... oraz przynaleznoscia wezlow do klastrow
					// z ta roznica, za w oryginalnej agregacji losowalismy do polaczenia etykiety (dostepne klastry), tutaj zas losujemy poszczegolne wezly,
					// ... co powoduje, ze mozliwe jest dokladanie polaczen w juz istniejacym klastrze, co z kolei implikuje, ze nie za kazdym razem nastapi
					// ... w ogole agregacja klastrow
					
					
					
					// 1. zeby sprawdzic, czy naleza do tego samego klastra, trzeba uzyskac informacje o przynaleznosci obu tych wezlow
					
					int przynaleznosc_wezla_i = membership_of_monomers.at(i);
					int przynaleznosc_wezla_j = membership_of_monomers.at(j);
					
					
					// 2. sprawdzic czy wysolosowane do polaczenia wezly i oraz j naleza do tego samego klastra - wtedy nic nie robimy,
					// ... jesli naleza do roznych klastrow, to trzeba zagregowac te klastry
					
					if( przynaleznosc_wezla_i != przynaleznosc_wezla_j )
					{
						// przynaleznosc klastra to to samo co etykieta klastra do zagregowania - przezwijmy to dla jasnosci
						
						int etykieta_klastra_1 = przynaleznosc_wezla_i;
						int etykieta_klastra_2 = przynaleznosc_wezla_j;
						
						// 3. dwa klastry o danych etykietach nalezy zagregowac
						
						// ind1 i ind2 to etykiety klastrow do zagregowania 
						// zagwarantujmy sobie zeby ind1 bylo zawsze mniejsze niz ind2
						if ( etykieta_klastra_2 < etykieta_klastra_1 )
						{
							int help = etykieta_klastra_1;
							etykieta_klastra_1 = etykieta_klastra_2;
							etykieta_klastra_2 = help;
						}

						int this_will_grow = etykieta_klastra_1;
						int this_will_vanish = etykieta_klastra_2;

						int ind_mer_1, ind_mer_2;

						// I have labels of mers, so here, I have to find the indexes of these mers in my memory vectors (labels and masses) to update masses
						// of course, labels and masses could be a tuple in one vector or sth
						for (unsigned int e = 0; e<labels_of_mers.size(); e++)
							if (labels_of_mers.at(e) == this_will_grow)
							    ind_mer_1 = e;
							else if (labels_of_mers.at(e) == this_will_vanish)
							    ind_mer_2 = e;

						 // In the vector of labels of mers, I delete the mer which is anihilated.
						 labels_of_mers.erase(labels_of_mers.begin() + ind_mer_2);

						 // I increase the mass of the mer which is growing
						 // and I delete the mass of the mer which is anihilated.
						 masses_of_mers.at(ind_mer_1) += masses_of_mers.at(ind_mer_2);
						 masses_of_mers.erase(masses_of_mers.begin() + ind_mer_2);

						 // Here, I change the membership of monomers. The monomers which were assigned to mer which vanished
						 // will be now assigned to the mer which is grown.
						 for (auto& a : membership_of_monomers)
						 	if (a == this_will_vanish)
							    a = this_will_grow;
					}
						    
						    
					// now, I have updated information on existing clusters (their labels - numbers), their masses and on the membership of the nodes
					// from the above I can later extract a giant cluster
					 
					// testowo wypiszmy sobie wielkosc labels of mers
					//cout << "labels of mers size: " << labels_of_mers.size() << endl;
					
					/**/
					
					
					
				}

				// odleglosc w przestrzeni ukrytej zapisywana jest zawsze
				//connectivity_matrix.at(j).at(i).second = dist;
				//connectivity_matrix.at(i).at(j).second = dist;

			}	
		} 



		 
		cout << "po generacji" << endl;
		cout << "licze rzeczywiste stopnie wierzcholkow" << endl;

		// sprawdzic czy rzeczywiscie sredni stopien jest <k> = 6

		for(int i=0; i<N; i++)
		{
		
			/*
			int k = 0;

			for(int j=0; j<N; j++)
			{
				// sumowac polaczenia danego wierzcholka
				//if ( connectivity_matrix.at(i).at(j).first == true )
				if ( check_connection( connectivity_list, i, j ) == true )
					k++;
			}

			nodes_.at(i).k = k;
			*/
			nodes_.at(i).k = get_degree(connectivity_list, i);
		}

		double srednie_k = 0.0;

		for(auto& n : nodes_)
			srednie_k += n.k;

		srednie_k /= double(N);

		cout << "srednie k rzeczywiste jest " << srednie_k << endl;
		cout << "komponentow rozlacznych: " << labels_of_mers.size() << endl;
		
		/*
		for(auto a:labels_of_mers)
			cout << a << " " ;
		cout << endl << endl;
		
		for(auto a:masses_of_mers)
			cout << a << " " ;
		cout << endl << endl;
		*/
		
		//for(auto a:membership_of_monomers)
		//	cout << a << " " ;
		//cout << endl << endl;

		int masa_najwiekszego_klastra = *max_element( masses_of_mers.begin(), masses_of_mers.end() );
		
		int indeks_najw_masy = 0;
		for(int i=0; i<masses_of_mers.size(); i++)
			if( masses_of_mers.at(i) > masses_of_mers.at(indeks_najw_masy) )
				indeks_najw_masy = i;
				
		int etykieta_najwiekszego_klastra = labels_of_mers.at( indeks_najw_masy );

		cout << "etykieta najw klastra " << etykieta_najwiekszego_klastra << endl;
		cout << "masa najw klastra " << masa_najwiekszego_klastra << endl;			

		

		// stworzmy osobna siec na przechowanie najwiekszego komponentu
		
		//int N_giant_comp = masa_najwiekszego_klastra;
		
		// w tej macierzy polaczyc bedziemy przechowywac tylo te polaczenia, ktore dotycza najwiekszego komponentu
		/*
		vector<vector<pair<bool,double>>> connectivity_matrix_giant_comp(N);

		for(int i=0; i<connectivity_matrix_giant_comp.size(); i++)
			connectivity_matrix_giant_comp.at(i) = vector<pair<bool,double>>(N, make_pair(false, 0.)); // we shall maybe use only the upper-right half

		for(int i=0; i<N; i++)
			for(int j=i+1; j<N; j++) 
			{	
				// sprawdzmy, czy wezly i oraz j uczestnicza w najwiekszym komponencie
				if( membership_of_monomers.at(i) == membership_of_monomers.at(j) == etykieta_najwiekszego_klastra )
				{
					// skoro tak, to sprawdzmy czy istnieje miedzy nimi polaczenie
					
					musimy przepisac je do nowej macierzy wylacznie dla najwiekszego klastra
					if( connectivity_matrix.at(i).at(j).first == true )
					{
						output_topology << i << " " << j << " " << connectivity_matrix.at(i).at(j).second << endl; 
					}
					
				}
				
			}
		/**/



		// histogram
		/*	
		vector<int> histogram(N);

		for(auto& n : nodes_)
		{
			histogram.at( n.k )++;
		}
		/**/

		//for(auto& h : histogram)
		//	cout << h << "  " ;
		//cout << endl;

		
		
		// zapiszmy do pliku wektor stopni wierzcholkow
		// oraz wektor stopni wierzcholkow dla najwiekszego klastra 
		/*
		string filename = "data_topologie/network_" + to_string(iteracje) + "_degrees.txt";
		ofstream output_degrees(filename.c_str());
		
		for(int i=0; i<nodes_.size(); i++)
		{
			output_degrees << nodes_.at(i).k << endl;
		}
		
		output_degrees.close();
		
		filename = "data_topologie/network_" + to_string(iteracje) + "_degrees_giantcluster.txt";
		ofstream output_degrees_gc(filename.c_str());
		
		for(int i=0; i<nodes_.size(); i++)
		{
			
			if( membership_of_monomers.at(i) == etykieta_najwiekszego_klastra )
				
				output_degrees_gc << nodes_.at(i).k << endl;
		}
		
		output_degrees_gc.close();
		/**/
		
		
		// zapiszmy do pliku wektor wierzcholkow
		
		cout << "zapisuje liste wierzcholkow z ich polozeniami na kole" << endl;
		
		
		string filename = where_to_save + "/network_" + to_string(iteracje) + "_nodes.txt";
		ofstream output_nodes(filename.c_str());
		
		//output_nodes << N << endl;
		//output_nodes << gamma << endl;
		//output_nodes << alfa << endl;
		
		output_nodes << setprecision(8);
		
		for(int i=0; i<nodes_.size(); i++)
		{
			//output_nodes << i << " " << nodes_.at(i).theta << " " << nodes_.at(i).kappa << " " << nodes_.at(i).k << endl;
			output_nodes << i << " " << nodes_.at(i).theta << endl;
		}
		
		output_nodes.close();
		/**/
		
		
		/*
		// zapiszmy macierz polaczen do pliku
		string filename_conmat = "network_" + to_string(iteracje) + "_connmatr.txt";
		ofstream output_connmatr(filename_conmat.c_str());
		
		for(int i=0; i<N; i++)
		{
			for(int j=0; j<N; j++)
			{
				output_connmatr << connectivity_matrix.at(i).at(j).first << " "; 
			}
			output_connmatr << endl;
		}
		
		output_connmatr.close();
		*/
		
		
		
		// zamiast powyzszych zapisow, zapiszmy siec w postaci "skad-dokad-dlugosc polaczenia"
		
		cout << "zapisuje liste polaczen z odlegosciami w przest. ukrytej H2" << endl;
		
		string filename_topology = where_to_save + "/network_" + to_string(iteracje) + "_topology_H2.edgelist";
		ofstream output_topology(filename_topology.c_str());
		
		//output_topology << N << endl;
		//output_topology << gamma << endl;
		//output_topology << alfa << endl;
		
		output_topology << setprecision(8);
		
		
		/*
		for(int i=0; i<N; i++)
		{
			for(int j=i+1; j<N; j++)
			{
			
				// sprawdzmy, czy wezly i oraz j uczestnicza w najwiekszym komponencie i tylko w tym wypadku je wypiszmy
				if( membership_of_monomers.at(i) == etykieta_najwiekszego_klastra and membership_of_monomers.at(j) == etykieta_najwiekszego_klastra )
				
				// wypisujemy do topologii wszystkie krawedzie obojetnie do jakiego klastra naleza, czyli wypisujemy do pliku cala siec !
				{
				
					//if( connectivity_matrix.at(i).at(j).first == true )
					//if( connectivity_matrix.at(i).at(j) == true )
					if( check_connection(connectivity_list, i, j) == true )
					{
						//output_topology << i << " " << j << " " << connectivity_matrix.at(i).at(j).second << endl;
						output_topology << i << " " << j << " " << d_theta_theta(nodes_.at(i).theta, nodes_.at(j).theta, R) << endl; 
						//cout << i << " " << j << " " << connectivity_matrix.at(i).at(j).second << endl; 
						
					}
				}
			}
		}
		*/
		
		
		
		// zmienne specjalnie dla tlumaczenia zmiennych ukrytych z S1 na H2
		
		// odnosnie kappa_0 zostawmy je na razie tak jak jest (po poprawce), ale nie wiem czy powinno ono byc brane przed poprawka, po poprawce
		// czy nie ma to znaczenia
		
		double R_daszek = 2*log(N/(mju*PI*kappa_0*kappa_0));
		
		
		
		// sposob zapisu sieci specjalnie dla listy polaczen
		
		unsigned node_i = 0;
		
		for(list<vector<unsigned>>::const_iterator it = connectivity_list.begin(); it != connectivity_list.end(); ++it)
		{
		
			// przeliczenie polozenia na przestrzen H2
			double radial_coord_i = R_daszek - 2*log( nodes_.at(node_i).kappa / kappa_0 );
		
		
			// sprawdzmy, czy wezel i uczestniczy w najwiekszym komponencie i tylko w tym wypadku wypiszmy polaczenia tego wezla
			// oczywiscie jesli mamy tryb zapisu calej sieci "fn", to bedziemy zapisywac wszystko
			if( what_to_save == "fn" or membership_of_monomers.at(node_i) == etykieta_najwiekszego_klastra )
			{	
				for(int j=0; j<(*it).size(); j++)
				{
					int node_j = (*it).at(j);
					//output_topology << node_i << " " << node_j << " " << d_theta_theta(nodes_.at(node_i).theta, nodes_.at(node_j).theta, R) << endl;
					//output_topology << node_i << " " << node_j << endl;
					
					
					// przeliczenie polozenia na przestrzen H2
					double radial_coord_j = R_daszek - 2*log( nodes_.at(node_j).kappa / kappa_0 );
					
					
					// obliczenie odleglosci hiperbolicznej w przestrzeni H2
					//double delta = d_theta_theta( nodes_.at(node_i).theta, nodes_.at(node_j).theta, R );
					double delta = delta_theta( nodes_.at(node_i).theta, nodes_.at(node_j).theta );
					//double odleglosc_hiperboliczna = radial_coord_i + radial_coord_j + 2*log( delta / 2. );
					double odleglosc_hiperboliczna = abs( radial_coord_i + radial_coord_j + 2*log( delta / 2. ) );
					/* // czasami (rzadko) zdarzaja sie odleglosci lekko < 0. 
					if(odleglosc_hiperboliczna<0.0)
					{
						cout << "odleglosc < 0, wynosi: " << odleglosc_hiperboliczna << endl;
						//getchar();
					}
					else if(odleglosc_hiperboliczna==0.0)
					{
						cout << "odleglosc == 0" << endl;
						getchar();
					}*/
					
					
					output_topology << node_i << " " << node_j << " " << odleglosc_hiperboliczna << endl;
				}
				
			}
			
			node_i++;
		}
		
		
		output_topology.close();
		/**/
		
		
		
		cout << "zapisuje liste polaczen z odlegosciami w przest. ukrytej S1" << endl;
		
		filename_topology = where_to_save + "/network_" + to_string(iteracje) + "_topology_S1.edgelist";
		output_topology = ofstream(filename_topology.c_str());
		
		output_topology << setprecision(8);
		
		
		// sposob zapisu sieci specjalnie dla listy polaczen
		
		node_i = 0;
		
		for(list<vector<unsigned>>::const_iterator it = connectivity_list.begin(); it != connectivity_list.end(); ++it)
		{
			// sprawdzmy, czy wezel i uczestniczy w najwiekszym komponencie i tylko w tym wypadku wypiszmy polaczenia tego wezla
			// oczywiscie jesli mamy tryb zapisu calej sieci "fn", to bedziemy zapisywac wszystko
			if( what_to_save == "fn" or membership_of_monomers.at(node_i) == etykieta_najwiekszego_klastra )
			{	
				for(int j=0; j<(*it).size(); j++)
				{
					int node_j = (*it).at(j);
					output_topology << node_i << " " << node_j << " " << d_theta_theta(nodes_.at(node_i).theta, nodes_.at(node_j).theta, R) << endl;
					//output_topology << node_i << " " << node_j << endl;
					
				}
				
			}
			
			node_i++;
		}
		
		
		output_topology.close();
		/**/
		
		
		// tutaj wykonamy specjalny zapis, tj. zapis wszystkich mozliwych polaczen w sieci z odleglosciami policzonymi wg S1 i wg H2
		// tzn de facto wszystkie odleglosci w przestrzeni ukrytej
		
		cout << "zapisuje liste wszystkich mozliwych polaczen z odlegosciami w przest. ukrytej S1" << endl;
		
		filename_topology = where_to_save + "/network_" + to_string(iteracje) + "_alldists_S1.edgelist";
		output_topology = ofstream(filename_topology.c_str());
		
		output_topology << setprecision(8);
		
		for(int node_i=0; node_i<N; node_i++)
		{
			if ( node_i%100 == 0 ) cout << " " << double(node_i)/N*100 << " % " << endl;
			
			for(int node_j=0; node_j<N; node_j++)
			{
				if(node_i == node_j) continue;
				
				output_topology << node_i << " " << node_j << " " << d_theta_theta(nodes_.at(node_i).theta, nodes_.at(node_j).theta, R) << endl;
			}
		}
		output_topology.close();
		
		
		
		
		cout << "zapisuje liste wszystkich mozliwych polaczen z odlegosciami w przest. ukrytej H2" << endl;
		
		filename_topology = where_to_save + "/network_" + to_string(iteracje) + "_alldists_H2.edgelist";
		output_topology = ofstream(filename_topology.c_str());
		
		output_topology << setprecision(8);
		
		for(int node_i=0; node_i<N; node_i++)
		{
			if ( node_i%100 == 0 ) cout << " " << double(node_i)/N*100 << " % " << endl;
		
			// przeliczenie polozenia na przestrzen H2
			double radial_coord_i = R_daszek - 2*log( nodes_.at(node_i).kappa / kappa_0 );
			
			for(int node_j=0; node_j<N; node_j++)
			{
				if(node_i == node_j) continue;
				
				// przeliczenie polozenia na przestrzen H2
				double radial_coord_j = R_daszek - 2*log( nodes_.at(node_j).kappa / kappa_0 );
				
				// obliczenie odleglosci hiperbolicznej w przestrzeni H2
				//double delta = d_theta_theta( nodes_.at(node_i).theta, nodes_.at(node_j).theta, R );
				double delta = delta_theta( nodes_.at(node_i).theta, nodes_.at(node_j).theta );
				//double odleglosc_hiperboliczna = radial_coord_i + radial_coord_j + 2*log( delta / 2. );
				double odleglosc_hiperboliczna = abs( radial_coord_i + radial_coord_j + 2*log( delta / 2. ) );
				
				output_topology << node_i << " " << node_j << " " << odleglosc_hiperboliczna << endl;
			}
		}
		output_topology.close();
		
		
		
		/*
		
		// odleglosci pomiedzy wszystkimi wezlami potrzebne są do greedy routing
		
		string filename_conmatrices = "data_conmatrices/network_" + to_string(iteracje) + "_conmatr.txt";
		ofstream output_conmatrices(filename_conmatrices.c_str());
		
		output_conmatrices << N << endl;
		output_conmatrices << gamma << endl;
		output_conmatrices << alfa << endl;
		
		output_conmatrices << setprecision(18);
		
		for(int i=0; i<N; i++)
		{
			for(int j=i+1; j<N; j++)
			{
				if( true )
				{
					output_conmatrices << i << " " << j << " " << connectivity_matrix.at(i).at(j).first << " " << connectivity_matrix.at(i).at(j).second << endl; 
					//cout << i << " " << j << " " << connectivity_matrix.at(i).at(j).second << endl; 
					
				}
			}
		}
		
		output_conmatrices.close();
		
		/*
		for(int i=0; i<N; i++)
		{
			for(int j=0; j<N; j++)
			{
				
					cout << connectivity_matrix.at(i).at(j).second << " ";
			}
			cout << endl;
		}
		/**/
		
		
		//dla kazdej sieci: dodatkowy plik zbiorczy jak zachowywalo sie ps dla tej realizacji sieci po uruchomieniu greedy routing
		
		//continue;


		// sprawdzic czy rezczywiscie rozklad k jest jak k^-gamma

		/*
		ofstream output_file("hist_out.txt");
		for(auto& h : histogram)
			output_file << h << endl;
		output_file.close(); 
		*/


		
		/*
		//int max_gr = 10e4; // w pracy 10e6
		int max_gr;
		//if(argc == 6)
		//	max_gr = atoi(argv[5]);
		//else
			max_gr = 10e4;


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
		
		/**/


		/*
		cout << "histo in bound" << endl;

		// histogramowanie stopni k wezlow
		
		vector<int> histogram(200);
		
		for(auto n : nodes_)
		{
			if(n.k<200)
				histogram.at(n.k)++;
		}
		
		// print histogram
		ofstream histo_file("histo.txt");
		for(auto h : histogram)
			histo_file << h << endl;
		histo_file.close();
		
		cout << "histo saved" << endl;
		
		// sprawdzilem empirycznie, ze liczyc slope gamma najbezipieczniej byloby od k=6 w gore gdyz dla wszystkich konfiguracji parametrow to zapewnia
		// ze liczymy slope tylko z tej czesci histogramu ktora na wykresie log-log rzeczywiscie jest prosta
		// ..ale jednak gammy wychodza za duze, modyfikuje wiec "prob odciecia" do k=4.
		// ..wydaje sie ze ten prob powinien byc dynamicznie dopasowywany w zaleznosci od z tego z jakim gamma generujemy (gdzie zaczyna sie slope?)
		
		// liczenie slope'a za pomoca rownania od MEJ Newman
		
		double gamma_z_danych = 0.0;
		double suma = 0.0;
		int liczba_probek = 0;
		
		for(auto n : nodes_)
			if(n.k >= 4)
			{
				suma += log( n.k / 4. );
				liczba_probek++;
			}
			
		gamma_z_danych = 1. + liczba_probek / suma;
		
		cout << "gamma z danych: " << gamma_z_danych << endl;
		
		getchar();
		
		/**/


	} ///// koniec for po instacjach sieci


	// obliczam srednie success ratio dla danego N, gamma i alfa
	/*
	double sredni_s_r = 0.0;

	for(auto& w:wyniki_gr_kontener_instacje)
		if(w > 0)
			sredni_s_r += w;

	sredni_s_r /= max_iteracji;

	cout << "sredni p_s dla N=" << N << ", gamma=" << gamma << ", i alfa=" << alfa << " jest " << sredni_s_r << endl;
	/**/
	
	

	


	//dla kazdego N, gamma, alfa : dodatkowy plik zbiorczy jak zachowywalo sie ps dla kolejnych realizacji sieci po uruchomieniu greedy routing
	/*
	string filename_grsummary = "grsummary_N=" + to_string_with_precision(N) + "_gamma=" + to_string_with_precision(gamma) + "_alfa=" + to_string_with_precision(alfa) + ".txt";
	ofstream output_grsummary(filename_grsummary.c_str());
		
	for(int i=0; i<max_iteracji; i++)
	{
		output_grsummary << i << " " << wyniki_gr_kontener_instacje.at(i) << endl; 
	
	}
		
	output_grsummary.close();
	/**/
	
	
	
	
	



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


	cout << "Wygenerowano " << max_iteracji << " sieci." << endl;





	return 0;
}
