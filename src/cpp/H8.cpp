/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "H8.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CH8::CH8()
{
    NEN_ = 8;	// Each element has 8 nodes
    nodes_ = new CNode*[NEN_];

    ND_ = 24;
    LocationMatrix_ = new unsigned int[ND_];

    ElementMaterial_ = nullptr;
}

//	Desconstructor
CH8::~CH8()
{
    if (nodes_)
        delete[] nodes_;

    if (LocationMatrix_)
        delete[] LocationMatrix_;
}

//	Read element data from stream Input
bool CH8::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int N;

    Input >> N;	// element number

    if (N != Ele + 1)
    {
        cerr << "*** Error *** Elements must be inputted in order !" << endl
             << "    Expected element : " << Ele + 1 << endl
             << "    Provided element : " << N << endl;

        return false;
    }

    unsigned int MSet;	// Material property set number
    unsigned int N1, N2,N3, N4, N5, N6, N7, N8;	// Left node number and right node number

    Input >> N1 >> N2 >> N3 >> N4 >> N5 >> N6 >> N7 >> N8 >> MSet;//8个节点（事实上存在nodes里）
    ElementMaterial_ = dynamic_cast<CH8Material*>(MaterialSets) + MSet - 1;
    nodes_[0] = &NodeList[N1 - 1];//存节点号的指针，注意nodes_在前面初始化过指向到数组，这里改为4个
    nodes_[1] = &NodeList[N2 - 1];
    nodes_[2] = &NodeList[N3 - 1];
    nodes_[3] = &NodeList[N4 - 1];
    nodes_[4] = &NodeList[N5 - 1];
    nodes_[5] = &NodeList[N6 - 1];
    nodes_[6] = &NodeList[N7 - 1];
    nodes_[7] = &NodeList[N8 - 1];

    for (int i = 0; i < 8; i++)
    {
        nodes_[i]->Tec_flag = true;
    }

	//!< Calculate nodal mass
	double rho = ElementMaterial_->density_0;
	//! two point Guass Quadrature
	double GP[2], weight[2], Jacobi, SHP[8];
	Guassian(2, GP, weight);
	double B[6][24]; // No use here

	for (int I = 0; I < 2; I++)
	{
		for (int J = 0; J < 2; J++)
		{
			for (int K = 0; K < 2; K++)
			{
				StrainMatrix(B, GP[I], GP[J], GP[K], Jacobi);
				double Gmass = rho * weight[I] * weight[J] * weight[K] * abs(Jacobi); // we can get volumn by integrating the function: f(x)=1
				// Expolate the mass from Gauss point to nodes
				SHPFunction(SHP, GP[I], GP[J], GP[K]);
				for (int i = 0; i < 8; i++)
				{
					nodes_[i]->mass += Gmass * SHP[i];
				}
			}
		}
	}
    return true;
}

//	Write element data to stream
void CH8::Write(COutputter& output, unsigned int Ele)
{
    output << setw(5) << Ele + 1 << setw(19) << nodes_[0]->NodeNumber
           << setw(9) << nodes_[1]->NodeNumber << setw(9) << nodes_[2]->NodeNumber << setw(9) << nodes_[3]->NodeNumber
           << setw(9) << nodes_[4]->NodeNumber << setw(9) << nodes_[5]->NodeNumber << setw(9) << nodes_[6]->NodeNumber
           << setw(9) << nodes_[7]->NodeNumber
           << setw(12) << ElementMaterial_->nset << endl;
}

void CH8::Write(CTecOutputter& output)
{
    output << nodes_[0]->NodeTecNumber << " " << nodes_[1]->NodeTecNumber << " " << nodes_[2]->NodeTecNumber << " " << nodes_[3]->NodeTecNumber << " "<< nodes_[4]->NodeTecNumber << " "
           << nodes_[5]->NodeTecNumber << " " << nodes_[6]->NodeTecNumber << " " << nodes_[7]->NodeTecNumber << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CH8::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; N++)
        for (unsigned int D = 0; D < 3; D++)
            LocationMatrix_[i++] = nodes_[N]->bcode[D];//LM，连接阵
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 8 node H8 element, element stiffness is a 24x24 matrix, whose upper triangular part
//	has 300 elements
unsigned int CH8::SizeOfStiffnessMatrix()
{
    return 300;
}

//	Calculate element stiffness matrix
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CH8::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());
    int Ngauss=2;
    double WGS[2]= {1.0,1.0},XGS[2]= {-0.577350269189625,0.577350269189625};
    double ks,yi,ph,wks,wyi,wph,detA;
    double Ke[24][24]= {0};
    //对所有高斯点循环
    for ( int i1=0; i1<Ngauss; i1++)
    {
        ks=XGS[i1];
        wks=WGS[i1];
        for ( int j1=0; j1<Ngauss; j1++)
        {
            yi=XGS[j1];
            wyi=WGS[j1];
            for ( int k1=0; k1<Ngauss; k1++)
            {
                ph=XGS[k1];
                wph=WGS[k1];
                //初始化并计算应变阵
                double B[6][24]= {0};
                StrainMatrix(B, ks,yi, ph, detA);



                //定义材料阵
                double Dmat[6][6]= {0};

                Constitutive(Dmat);

                //计算应力阵S=DBT
                double S[6][24]= {0};
                for(unsigned int N1=0; N1<24; N1++)
                {
                    for(unsigned int N2=0; N2<6; N2++)
                    {
                        for(unsigned int N3=0; N3<6; N3++)
                        {
                            S[N2][N1]+=Dmat[N2][N3]*B[N3][N1];
                        }
                    }
                }

                //计算通过高斯积分计算刚度阵
                double WT=abs(detA)*WGS[i1]*WGS[j1]*WGS[k1];
                for(unsigned int N1=0; N1<24; N1++)
                {
                    for(unsigned int N2=0; N2<24; N2++)
                    {
                        for(unsigned int N3=0; N3<6; N3++)
                        {
                            Ke[N1][N2]+=B[N3][N1]*S[N3][N2]*WT;
                        }
                    }
                }

                int N3=0;
                for (int N1=0; N1<24; N1++)
                {
                    for(int N2=N1; N2>=0; N2--)
                    {
                        Matrix[N3]=Ke[N2][N1];
                        N3++;
                    }
                }

            }
        }
    }
}


//*****************************************************************************************************************/
//	Calculate element stress

void CH8::ElementStress(double* stress, double* Displacement)
{
    CH8Material* material_ = dynamic_cast<CH8Material*>(ElementMaterial_);	// Pointer to material of the element



    int Ngauss=2;
    double WGS[2]= {1.0,1.0},XGS[2]= {-0.577350269189625,0.577350269189625};
    double ks,yi,ph,wks,wyi,wph,detA;
    double Ke[24][24]= {0};
    //先将应力回0
    clear(stress, 48);

    for ( int i1=0; i1<Ngauss; i1++)
    {
        ks=XGS[i1];
        wks=WGS[i1];
        for ( int j1=0; j1<Ngauss; j1++)
        {
            yi=XGS[j1];
            wyi=WGS[j1];
            for ( int k1=0; k1<Ngauss; k1++)
            {
                ph=XGS[k1];
                wph=WGS[k1];
                //初始化并计算应变阵
                double B[6][24]= {0};
                StrainMatrix(B, ks,yi, ph, detA);



                //定义材料阵
                double Dmat[6][6]= {0};

                Constitutive(Dmat);

                //计算应力阵S=DBT
                double S[6][24]= {0};
                for(unsigned int N1=0; N1<24; N1++)
                {
                    for(unsigned int N2=0; N2<6; N2++)
                    {
                        for(unsigned int N3=0; N3<6; N3++)
                        {
                            S[N2][N1]+=Dmat[N2][N3]*B[N3][N1];
                        }
                    }
                }

//***********************************************************************
                //计算应力stress=S*U
                for (int k = 0; k < 24; k++)
                {
                    if (LocationMatrix_[k])
                    {
                        stress[((i1 * 2 + j1) * 2 + k1) * 6] += S[0][k] * Displacement[LocationMatrix_[k] - 1];
                        stress[((i1 * 2 + j1) * 2 + k1) * 6 + 1] += S[1][k] * Displacement[LocationMatrix_[k] - 1];
                        stress[((i1 * 2 + j1) * 2 + k1) * 6 + 2] += S[2][k] * Displacement[LocationMatrix_[k] - 1];
                        stress[((i1 * 2 + j1) * 2 + k1) * 6 + 3] += S[3][k] * Displacement[LocationMatrix_[k] - 1];
                        stress[((i1 * 2 + j1) * 2 + k1) * 6 + 4] += S[4][k] * Displacement[LocationMatrix_[k] - 1];
                        stress[((i1 * 2 + j1) * 2 + k1) * 6 + 5] += S[5][k] * Displacement[LocationMatrix_[k] - 1];
                    }
                }
            }
        }
    }

    double str[6];
    str[0] = stress[0] + stress[6] + stress[12] + stress[18] + stress[24] + stress[30] + stress[36] + stress[42];
    str[1] = stress[1] + stress[7] + stress[13] + stress[19] + stress[25] + stress[31] + stress[37] + stress[43];
    str[2] = stress[2] + stress[8] + stress[14] + stress[20] + stress[26] + stress[32] + stress[38] + stress[44];
    str[3] = stress[3] + stress[9] + stress[15] + stress[21] + stress[27] + stress[33] + stress[39] + stress[45];
    str[4] = stress[4] + stress[10] + stress[16] + stress[22] + stress[28] + stress[34] + stress[40] + stress[46];
    str[5] = stress[5] + stress[11] + stress[17] + stress[23] + stress[29] + stress[35] + stress[41] + stress[47];



    for (unsigned int N = 0; N < NEN_; N++)
    {
        nodes_[N]->count += 8;
        nodes_[N]->Stress[0] += str[0];		// stress xx
        nodes_[N]->Stress[1] += str[1];		// stress yy
        nodes_[N]->Stress[2] += str[2];		// stress zz
        nodes_[N]->Stress[3] += str[3];		// stress xy
        nodes_[N]->Stress[4] += str[4];		// stress xz
        nodes_[N]->Stress[5] += str[5];		// stress yz
    }

}

void CH8::Constitutive(double Dmat[6][6])
{


    CH8Material* material_ = dynamic_cast<CH8Material*>(ElementMaterial_);	// Pointer to material of the element
    double EMoudule=material_->E,poi=material_->Poisson;
    double namda=EMoudule*poi/(1+poi)/(1-2*poi),GMoudule=EMoudule/2/(1+poi);
    for (unsigned int N1=0; N1<3; N1++)
    {
        for(unsigned int N2=0; N2<3; N2++)
        {
            if(N1==N2)
            {
                Dmat[N1][N2]=namda+2*GMoudule;
            }
            else
            {
                Dmat[N1][N2]=namda;
            }
        }
    }

    for (unsigned int N1=3; N1<6; N1++)
    {
        Dmat[N1][N1]=2*GMoudule;
    }
}

//! Return the shape function value of point with parent coordinate (xi, eta)
void CH8::SHPFunction(double SHP[8], double R, double S , double T)
{
	SHP[0]=(1-R)*(1-S)*(1-T);
	SHP[1]=(1+R)*(1-S)*(1-T);
	SHP[2]=(1+R)*(1+S)*(1-T);
	SHP[3]=(1-R)*(1+S)*(1-T);
	SHP[4]=(1-R)*(1-S)*(1+T);
	SHP[5]=(1+R)*(1-S)*(1+T);
	SHP[6]=(1+R)*(1+S)*(1+T);
	SHP[7]=(1-R)*(1+S)*(1+T);
}

void  CH8::StrainMatrix(double B[6][24], double ks, double yi, double ph,double &detA)
{
    /*按(-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),
    (-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1)的顺序*/
    //设循环体也可以
    double dndk[8][3];
    dndk[0][0]=-(1-yi)*(1-ph)/8;
    dndk[0][1]=-(1-ks)*(1-ph)/8;
    dndk[0][2]=-(1-yi)*(1-ks)/8;

    dndk[1][0]=(1-yi)*(1-ph)/8;
    dndk[1][1]=-(1+ks)*(1-ph)/8;
    dndk[1][2]=-(1-yi)*(1+ks)/8;

    dndk[2][0]=(1+yi)*(1-ph)/8;
    dndk[2][1]=(1+ks)*(1-ph)/8;
    dndk[2][2]=-(1+yi)*(1+ks)/8;

    dndk[3][0]=-(1+yi)*(1-ph)/8;
    dndk[3][1]=(1-ks)*(1-ph)/8;
    dndk[3][2]=-(1+yi)*(1-ks)/8;

    dndk[4][0]=-(1-yi)*(1+ph)/8;
    dndk[4][1]=-(1-ks)*(1+ph)/8;
    dndk[4][2]=(1-yi)*(1-ks)/8;

    dndk[5][0]=(1-yi)*(1+ph)/8;
    dndk[5][1]=-(1+ks)*(1+ph)/8;
    dndk[5][2]=(1-yi)*(1+ks)/8;

    dndk[6][0]=(1+yi)*(1+ph)/8;
    dndk[6][1]=(1+ks)*(1+ph)/8;
    dndk[6][2]=(1+yi)*(1+ks)/8;

    dndk[7][0]=-(1+yi)*(1+ph)/8;
    dndk[7][1]=(1-ks)*(1+ph)/8;
    dndk[7][2]=(1+yi)*(1-ks)/8;



    //初始化
    double dxdk[3][3]= {0},dkdx[3][3]= {0};


    //矩阵乘法求dxdk
    for( int N1=0; N1<3; N1++)
    {
        for( int N2=0; N2<3; N2++)
        {
            for( int N3=0; N3<8; N3++)
            {
                dxdk[N1][N2]+=dndk[N3][N2]*(nodes_[N3]->XYZ[N1]);
            }
        }
    }
    //下面对dxdk求逆
    detA=dxdk[0][0]*(dxdk[1][1]*dxdk[2][2]-dxdk[1][2]*dxdk[2][1])-dxdk[0][1]*(dxdk[1][0]*dxdk[2][2]-dxdk[1][2]*dxdk[2][0])+dxdk[0][2]*(dxdk[1][0]*dxdk[2][1]-dxdk[1][1]*dxdk[2][0]);
    
	if (detA == 0)
	{
		cout << "The Jacobi of an element equals zero!" << endl;
		exit(1);
	}
	
	
	int k11,k12,k21,k22;
    for( int N1=0; N1<3; N1++)
    {
        for( int N2=0; N2<3; N2++)
        {
            switch(N1)
            {
                case 0:
                    k11=1,k12=2;
                    break;
                case 1:
                    k11=0,k12=2;
                    break;
                case 2:
                    k11=0,k12=1;

            }
            switch(N2)
            {
                case 0:
                    k21=1,k22=2;
                    break;
                case 1:
                    k21=0,k22=2;
                    break;
                case 2:
                    k21=0,k22=1;

            }
            dkdx[N2][N1]=(dxdk[k11][k21]*dxdk[k12][k22]-dxdk[k11][k22]*dxdk[k12][k21])*pow(-1,(N1+N2+2))/detA;

        }
    }

    double dndx[8][3]= {0};
    for(unsigned int N1=0; N1<8; N1++)
    {
        for(unsigned int N2=0; N2<3; N2++)
        {
            for(unsigned int N3=0; N3<3; N3++)
            {
                dndx[N1][N2]+=dndk[N1][N3]*dkdx[N3][N2];
            }
        }
    }
    //生成B

    for (unsigned int N1=0; N1<8; N1++)
    {
        B[0][3*N1]=dndx[N1][0];
        B[1][3*N1+1]=dndx[N1][1];
        B[2][3*N1+2]=dndx[N1][2];

        B[3][3*N1]=dndx[N1][1];
        B[3][3*N1+1]=dndx[N1][0];

        B[4][3*N1]=dndx[N1][2];
        B[4][3*N1+2]=dndx[N1][0];

        B[5][3*N1+1]=dndx[N1][2];
        B[5][3*N1+2]=dndx[N1][1];
    }

}