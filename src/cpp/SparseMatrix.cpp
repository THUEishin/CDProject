#include "../h/SparseMatrix.h"
#include "../h/Domain.h"
#include <vector>

using namespace std;

//! constructors
template <class T>
CSparseMatrix<T>::CSparseMatrix()
{
	CDomain* FEMData = CDomain::Instance();
	unsigned int NUMNP = FEMData->GetNUMNP();
	unsigned int NUMEG = FEMData->GetNUMEG();
	CNode* nodelist = FEMData->GetNodeList();
	NEQ_ = FEMData->GetNEQ();

	vector<int> node_flag;
	vector<int> node_c_node;
	vector<int> node_c_freedom;
	vector<vector<int>> freedom_c_freedom;

	_iK = new unsigned int[NEQ_ + 1];
	_iK[0] = 1;
	freedom_c_freedom.resize(NEQ_);

//! Initialize the node flag
	node_flag.resize(NUMNP);
	for (unsigned int N = 0; N < NUMNP; N++)
	{
		node_flag[N] = 0;
	}

//! Calculate the cells connected to a node
	vector<vector<CElement*>> NE;
	NE.resize(NUMNP);
	for (unsigned int NEG = 0; NEG < NUMEG; NEG++)
	{
		CElementGroup& Elegrp = FEMData->GetEleGrpList()[NEG];
		unsigned int NUME = Elegrp.GetNUME();
		for (unsigned int E = 0; E < NUME; E++)
		{
			unsigned int NEN = Elegrp[E].GetNEN();
			CNode** nodes = Elegrp[E].GetNodes();
			for (unsigned int N = 0; N < NEN; N++)
			{
				NE[nodes[N]->NodeNumber - 1].push_back(&Elegrp[E]);
			}
		}
	}

//! Calculate freedom_c_freedom
	for (unsigned int N; N < NUMNP; N++)
	{
		unsigned int NDF = CNode::NDF;
		unsigned int flag= 0;
		//! Identify if the node has freedom
		for (unsigned int I = 0; I < NDF; I++)
		{
			flag += nodelist[N].bcode[I];
		}

		if (flag)
		{
			//If there is a freedom at the node, calculate the freedom id connected to this node

			//Calculate the node_c_node for node i
			int num_connected_elem = NE[N].size();
			for (int E = 0; E < num_connected_elem; E++)
			{
				CElement* Ele = NE[N][E];
				unsigned int NEN = Ele->GetNEN();
				CNode** nodes = Ele->GetNodes();
				for (unsigned int nen = 0; nen < NEN; nen++)
				{
					if (node_flag[nodes[nen]->NodeNumber - 1] == 0)
					{
						node_c_node.push_back(nodes[nen]->NodeNumber);
						node_flag[nodes[nen]->NodeNumber - 1] = 1;
					}
				}
			}

			//Clear node_flag after obtain the node_c_node
			for (int E = 0; E < num_connected_elem; E++)
			{
				CElement* Ele = NE[N][E];
				unsigned int NEN = Ele->GetNEN();
				CNode** nodes = Ele->GetNodes();
				for (unsigned int nen = 0; nen < NEN; nen++)
				{
					node_flag[nodes[nen]->NodeNumber - 1] = 0;
				}
			}

			//Order the node_c_node
			sort(node_c_node.begin(), node_c_node.end());

			//Calulate the id of freedom degree connected to this node
			int num_connected_node = node_c_node.size();
			for (int N = 0; N < num_connected_node; N++)
			{
				int node_id = node_c_node[j] - 1;
				for (unsigned int k = 0; k < NDF; k++)
				{
					if (nodelist[node_id].bcode[k] > 0)
						node_c_freedom.push_back(nodelist[node_id].bcode[k]);
				}
			}

			//Calculate the freedom_c_freedom of the freedom on this node
			for (unsigned int j = 0; j < NDF; j++)
			{
				int num_connected_freedom = node_c_freedom.size();
				if (nodelist[N].bcode[j] > 0)
				{
					for (unsigned int k = 0; k < num_connected_freedom; k++)
					{
						if (node_c_freedom[k] >= nodelist[N].bcode[j])
							freedom_c_freedom[nodelist[N].bcode[j] - 1].push_back(node_c_freedom[k]);
					}
				}
			}
		}
		node_c_node.clear();
		node_c_freedom.clear();
	}

	//Calcualte _iK and number of non-zero element
	int non_zero = 0;
	for (unsigned int i = 0; i < NEQ_; i++)
	{
		_iK[i + 1] = _iK[i] + freedom_c_freedom[i].size();
		non_zero += freedom_c_freedom[i].size();
	}

	//Calculte _jK
	_jK = new int[non_zero];
	non_zero = 0;
	for (unsigned int i = 0; i < NEQ_; i++)
	{
		int num_freedom_connected_freedom = freedom_c_freedom[i].size();
		for (int j = 0; j < num_freedom_connected_freedom; j++)
		{
			_jK[non_zero + j] = freedom_c_freedom[i][j];
		}
		non_zero = non_zero + num_freedom_connected_freedom;
	}
	NNZ_ = non_zero;

	//Allocate _K
	_K = new T[NNZ_];
	for (int i = 0; i < NNZ_; i++)
		_K[i] = T(0);
}

//! deconstructor
template <class T>
CSparseMatrix<T>::~CSparseMatrix()
{
	delete[] _K;
	delete[] _iK;
	delete[] _jK;
}