// C++ implementation to find the cantor set
// for n levels and
// for a given start_num and end_num

#include <bits/stdc++.h>
#include <vector>



//using namespace std;

// The Linked List Structure for the Cantor Set
typedef struct cantor {
	double start, end;
	struct cantor* next;
} Cantor;

// Function to initialize the Cantor Set List
Cantor* startList(Cantor* head,
				double start_num,
				double end_num)
{
	if (head == NULL) {
		head = new Cantor;
		head->start = start_num;
		head->end = end_num;
		head->next = NULL;
	}
	return head;
}

// Function to propagate the list
// by adding new nodes for the next levels
Cantor* propagate(Cantor* head)
{
	Cantor* temp = head;

	if (temp != NULL) {
		Cantor* newNode
			= new Cantor;
		double diff
			= (((temp->end) - (temp->start)) / 3);

		// Modifying the start and end values
		// for the next level
		newNode->end = temp->end;
		temp->end = ((temp->start) + diff);
		newNode->start = (newNode->end) - diff;

		// Changing the pointers
		// to the next node
		newNode->next = temp->next;
		temp->next = newNode;

		// Recursively call the function
		// to generate the Cantor Set
		// for the entire level
		propagate(temp->next->next);
	}

	return head;
}

// Function to print a level of the Set
void print(Cantor* temp)
{
	while (temp != NULL) {
		printf("[%lf] -- [%lf]\t",
			temp->start, temp->end);
		temp = temp->next;
	}
	std::cout << std::endl;
}

// Function to build and display
// the Cantor Set for each level
std::vector<double> buildCantorSetOnCircle(int A, int B, int L)
{
	std::vector<double> wektor;

	Cantor* head = NULL;
	head = startList(head, A, B);
	for (int i = 0; i < L; i++) {
		//cout <<"Level_"<< i<<" : ";
		//print(head);
		propagate(head);
	}
	//cout <<"Level_"<< L<<" : ";
	//print(head);
	
	auto temp = head;
	
	while (temp != NULL)
	{
		wektor.push_back(temp->start);
		wektor.push_back(temp->end);
		temp = temp->next;
	}
	
	double d_pi = 2 * 3.14159265358979323846;
	
	// transformacja na kolo
	for(auto& w:wektor)
		w *= d_pi;
	
	
	return wektor;
}





/*
// Driver code
int main()
{
	int A = 0;
	int B = 1;
	int L = 9;
	vector<double> zbior_cantora = buildCantorSetOnCircle(A, B, L);

	for(auto e:zbior_cantora)
		cout << e << endl;
		
	cout << "liczba odcinkow:" << zbior_cantora.size()/2 << endl;

	return 0;
}

// This code is contributed by shivanisingh
*/
