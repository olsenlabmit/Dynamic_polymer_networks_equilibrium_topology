// KMC program
/*
#######################################
#                                     #
#-- KMC for Dynamic Networks --#
#------  Author: Devosmita Sen  --------#
#                                     #
#######################################
// Updated Devosmita Sen Februrary 16, 2023
// KMC for dynamic networks
// A4-B2 code--- A2B4 used (Bifunctional polymer chain and tetrafunctional junctions)
//Step 1: only forward reactions (same as earlier code by Tzyy-Shyang Lin)
//Step 2: both forward abd reverse reactions (depends on probability of forward and reverse rxns)
//Step 3: Only forward reaction (same as Step 1, but starting at the point where Step 2 ends)

// Overview of reverse reaction implementation: 

select junction - JA and JB- based on random scission criterion- select a bond randomly and cleave it
update urA, urB, reactedA and reactedB arrays- these arrays contain the list of A and B junctions which are completely reacted
(in reactedA and reactedB) and those which have at least one unreacted site
sort urA and urB arrays- this is needed because this array is used later- so best to sort it wrt junction numbers
update loop- updates the loop fraction
update neighA and neighB- these are the arrays containing the neighbours of each A and B junction- 
update junctiondist- update the junction distance between each A and B junctions- 
this step is the main bottleneck- because once a bond is broken- all info about its shortest dist is lost, and hence
all distances have to be recalculated
the algo which is used in Floyd-Warshall algorithm- but this is applied only to molecules within the cluster 
because for all other molecules (outside that cluster to which JA an JB belonged to)- the distance between junctions 
will stay as infinity (used as 0 in the code!- as a conventiosn)
collect connected-AcA,AcB, BcA, BcB- update the A, B junctions connected to JA and JB respectively
update molecule number- update the cluster numbers after running
update sum- update sum according to updated junction distances
*/
#include<iostream>
#include<cstdio>
#include<vector>
#include<algorithm>
#include<cmath>
#include<ctime>
#include<omp.h>
#include<map>
#include<utility>
#include<queue>
#include<algorithm>
#include<fstream>
#include<string>
#include <numeric>
#include<stdexcept>


#define MATRIX_A 0x9908b0dfUL //for RNG
#define UPPER_MASK 0x80000000UL
#define LOWER_MASK 0x7fffffffUL
#define DIST_TYPE short
using namespace std;
// Global variables
/*********************************/
/* Simulation Parameters */
/*********************************/
// chemistry related parameters


double R2_func()
{
  double n;
  ifstream infile("R2.txt");
  if (infile>>n)
    return n;
  else
    throw std::runtime_error("Cannot read size from file");
}
const double Nbsq = R2_func();


size_t NRA_func()
{
  size_t n;
  ifstream infile("NRA.txt");
  if (infile>>n)
    return n;
  else
    throw std::runtime_error("Cannot read size from file");
}
const size_t NRA = NRA_func();





double c0_in_func() // input C0
{
  double n;
  ifstream infile("C0.txt");
  if (infile>>n)
    return n;
  else
    throw std::runtime_error("Cannot read size from file");
}
const double c0_in = c0_in_func();

const double ca0=c0_in*2*0.001; // can be overwritten w/ argument input into main

double invKeq_factor_func()
{
  double n;
  ifstream infile("invKeq.txt");
  if (infile>>n)
    return n;
  else
    throw std::runtime_error("Cannot read size from file");
}
const double invKeq_factor = invKeq_factor_func();
double invKeq=invKeq_factor; //*ca0;

const size_t wrt_step=NRA/10;


// this definition of conc is number of functional groups / volume- 0.780- is from paper- and the factor 2 helps to convert into KMC definition
// 0.001 factor to convert from mM to M
const double c_star = ca0*0.5 * (6.02214e23/1.0e24)*pow(Nbsq,1.5);
const size_t MAX=(2<<sizeof(DIST_TYPE)*8-1); // maximal loop size considered
// /* A4B2 settings (comment this line to use this setting)
const size_t fA=2,fB=4; //changed to A2B2 now
const char dA=1,dB=0;
//const size_t molsizeA=1,molsizeB=0;
size_t step;

class bond {       // The class
  public:             // Access specifier
    size_t JA;        // Junction A
    size_t JB;  //Junction B
};

// */
/* A4B4 settings (comment this line to use this setting)
const size_t fA=4,fB=4;
const char dA=1,dB=1;
const size_t mA=1,mB=1;
// */
// basic simulation settings (final conversion, system size, loop?, etc.)
double conversion=1.0; // final conversion
//size_t NRA=7500; // number of A precursor
// can be overwritten w/ argument input into main
const size_t NRB=NRA*fA/fB; // number of B precursor // changes dewpending on NA
const size_t NA=NRA*fA,NB=NRB*fB;
size_t start = NA>NB ? NB : NA;
size_t i_rxn=start;// index for tracking conversion
const bool LOOPFREE=false; // loop free or not, affects probability and final conversion
bool DEBUG = false; // debug mode, if TRUE, RNG will use the TEST_SEED specified
size_t TEST_SEED = 100; //
// output/environment settings
char PATH[]="";
size_t suffix;
bool writeLoop = true;
bool writeMW = false;
bool writeDistAA = false;
bool logDPdist = false;
bool logNW = false;
bool write_conv_loop=true;
bool write_connectivity_data=true;
char prefixLoopFrac[]="lpfrac_", prefixMW[]="MW_", prefixMWdist[]="cumDist_", prefixDegree[]="deg_", prefixDist[]="Dist_", prefixconv_loop_step1[]="1_Conv_loop_",prefixconv_loop_step2[]="2_Conv_loop_",prefixconv_loop_step3[]="3_Conv_loop_",prefixconnect_1[]="1_network_",prefixconnect_2[]="2_network_",prefixconnect_3[]="3_network_";
// logging parameters
size_t MWfreq=1; // output/record frequency for DPw-ie. output after every 10 steps
size_t loopVecSize = 2; // max order of loop to log
//size_t nNW = 2; // number of times to output network topology// why 2 times?
//size_t nDPdist = 2; // number of times to log MW distributions // why 2 times?
// RNG related
unsigned long mt[624]; // someything to do with RNG??
int mti=625;
double fn[128], wn[128];
int kn[128];
void RNG_initialize(unsigned long);
unsigned long rand32();
size_t seed; 
/********************************/
/* Simulation Variables */
/********************************/

//double sumrA_check=fA*NRA+1;
// variables relating to KMC
double PAB,la; // KMC propensities

size_t tmpJA,tmpJB; // temporary containers for junction numbers
vector<double> p; // probabilities of different reactions

//vector<size_t> G[NRB]; // adjacency list wrt only crosslinkers- required for shortest distance calculation in BFS
vector<vector<size_t> > G;
vector<vector<unsigned short> > JunctionDist; // array for pairwise topological distance
vector<vector<unsigned short> > LoopSize ; // array for pairwise loopsize-> if LoopSize[JA][JB]=1-> this means that JA and JB are involved in formation of a loop of size 1
vector<double> Conv; // conversions at which the following values were logged
vector<double> time_arr; 
vector<double> full_reactedA_array; // number of A junctions fully reacted
vector<double> all_reactedA_fg; // number of all A functional groups (not junctions) which are reacted
vector<double> Conv_calc; // conversions calculated using formula: conv=1-urA.size()/NA
vector<double> Sum; // stores the value of sum at each step // for dissoc

vector<bond> all_bonds;

vector<vector<size_t> > neighA,neighB; // neighbors of A(B) precursors

vector<size_t> urA,urB; // unreacted A and B junctions, size changing
vector<size_t> reactedA(0);//,reactedB; // unreacted A and B junctions, size changing
vector<size_t> reactedB(0);

vector<double> sumA; // sum of probability relating to each A junction// dissoc


double sum; // sum of propensities// dissociative

// variable related to rejection condition
double V,BoxSize; // V = volume of simulation box (nm^3)
// BoxSize = V^(1/3) in unit of nm
//vector<vector<double> > rA,rB; // spatial coordinates of A and B molecules

vector<size_t> mol; // molecule idx for each junction
// size = number of Junctions (NRA+NRB)
size_t largestMol; // index of largest molecule
clock_t c0;

//vector<double> loop0; //
vector<double> loop1; //
vector<double> loop2; //
vector<double> loop3; //
vector<double> loop4; //
/*vector<double> loop5; //
vector<double> loop6; //
vector<double> loop7; //
vector<double> loop8; //
vector<double> loop9; //
vector<double> loop10; //
*/
// Connectivity data for network

vector<size_t> node1; 
vector<size_t> node2; 
// all with size = NRA*fA reserved
// containers for logging number of loops for loop counting at end of KMC
vector<double> loop; // number of loops- if there is fA=2 or fB=2, then this will be 2*num_loops
vector<double> loopfrac; // loop fraction


/*********************/
/* Functions */
/*********************/
// KMC functions
void initialize();
bool KMCstep_forward(); // implements KMC step- forward
bool KMCstep_reverse(); // implements KMC step- reverse

double KMCconv(double,double ,double , double,size_t ); // outputs the conversion- i think this conversion is same as the actual conversion
void output(size_t);
void output_time_step(size_t, size_t);
// helper functions
// general helper functions
size_t dist(unsigned char);
int getJunctionDistAA(size_t,size_t); //part of the shortest path update step??
int getJunctionDistBB(size_t,size_t);
// steps in KMCstep()
void SelectJunct_forward(size_t &,size_t &,size_t &,size_t &);// pair selection algorithm 
void SelectJunct_reverse(size_t &,size_t &);//,size_t &,size_t &);// pair selection algorithm 

void UpdateConnectivityData_forward(size_t &JA,size_t &JB); // update connectivity data (B4 node 1 connected to B4 node 2)
void UpdateConnectivityData_reverse(size_t &JA,size_t &JB); // update connectivity data (B4 node 1 connected to B4 node 2)


void UpdateSum_forward(const size_t,const size_t,const size_t,const size_t); // propensity sum update

void UpdateSum_reverse(const size_t,const size_t);//,const size_t,const size_t); // propensity sum update
void UpdateLoop_forward(const size_t,const size_t);// loop count update
void UpdateLoop_reverse(const size_t,const size_t);// loop count update

void CollectConnected(const size_t,const size_t,vector<size_t>&,vector<size_t>&,vector<size_t>&,vector<size_t>&,vector<size_t>&,vector<size_t>&);
void CollectConnected_reverse(const size_t,const size_t,vector<size_t>&,vector<size_t>&,vector<size_t>&,vector<size_t>&,vector<size_t>&,vector<size_t>&);

void UpdateJuncDist_forward(const size_t,const size_t,const vector<size_t>&,const vector<size_t>&);
void UpdateJuncDist_reverse(const size_t,const size_t);//,const vector<size_t>&,const vector<size_t>&,const vector<vector<size_t> >&,const vector<vector<size_t> >&);
//void BFS_on_node(vector<size_t>  , vector<vector<unsigned short> > &, vector<vector< int> > &, vector <size_t> &);



void UpdateJuncDist_common(const size_t,const size_t,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&);//part of the shortest path update step??
//void UpdateJuncDist_reverse(const size_t,const size_t,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&);//part of the shortest path update step??
void UpdateJuncDist_common_reverse(const size_t,const size_t,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&);//part of the shortest path update step??

void UpdateMol_forward(const size_t,const size_t,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&);

void UpdateMol_reverse(const size_t,const size_t,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&);

// used in KMCconv()
void updateWriteData(double,double,size_t, double, double);

// MAIN FUNCTION
int main(int argc,char* argv[])// possible input arguements are- cstar, NRA, suffix
{
  //freopen("log_file", "w", stdout ); // opens file named "filename" for output 
    cout<<"invKeq"<<invKeq<<endl;
    //system("pause");

   
 
    // handle argument input
    /*if(argc<2) { }
    else {
    double c_star0 = c_star;
    c_star = atof(argv[1]);// atof- string to float
    ca0 = ca0 * (c_star/c_star0);
    }
    if(argc==3) {
    NRA = atoi(argv[2]);// atoi- string to integer
    NRB = atoi(argv[2])*fA/fB;
    } 
    else if(argc==4) {
    NRA = atoi(argv[2]);
    NRB = atoi(argv[2])*fA/fB;
    suffix = atoi(argv[3]); // this is a global variable, so it is not being used in main(), but still can be used in initialize() function 
    }
    */
    // log run time
    clock_t c;
    c=clock();
    //c0=clock();
    // RUN KMC simulation
    initialize();
    //cout<<"initialized\n";
    //double invKeq=1.0;
    
    //step 1
   // cout<<"step1\n";
    //cout<<"urA.size()"<<urA.size();
   // cout<<"reactedA.size()"<<reactedA.size();

    step=1;
    Conv.clear();
    time_arr.clear();
    full_reactedA_array.clear();
    all_reactedA_fg.clear();
    Conv_calc.clear();
    Sum.clear();
    
    //loop0.clear();
    loop1.clear();
    loop2.clear();
    loop3.clear();
    loop4.clear();
    /*
    loop5.clear();
    loop6.clear();
    loop7.clear();
    loop8.clear();
    loop9.clear();
    loop10.clear();*/
    
    node1.clear();
    node2.clear();

    double K_forward=1;
    double finalConv1 = KMCconv(conversion,0,K_forward,0,step);
    cout<<"KMCconv done\n";
    output(step);
    cout<<"Hi\n";
    //FILE * fp;
    //fp=fopen("molecule_size1.csv","a");

    cout<<"c_star\n";
    printf("%.8f\t",c_star);
    cout<<"\n";
    cout<<"loop_frac\n";
    for(size_t l=0;l<loop.size();++l) printf("%.8f\t",loopfrac[l]);
    cout<<"\n";
    cout <<"finalConv1 "<< finalConv1 << "\n";

    cout<< "BoxSize " << BoxSize << "\n";

    cout << "seed " <<seed << "\n";
    cout<< (double)(clock()-c)/CLOCKS_PER_SEC<<endl;


    //system("pause");
    //step 2
    cout<<"urA.size()"<<urA.size();
    cout<<"reactedA.size()"<<reactedA.size();

    Conv.clear();
    time_arr.clear();
    full_reactedA_array.clear();
    all_reactedA_fg.clear();
    Conv_calc.clear();
    Sum.clear();
    
    //loop0.clear();
    loop1.clear();
    loop2.clear();
    loop3.clear();
    loop4.clear();
    /*
    loop5.clear();
    loop6.clear();
    loop7.clear();
    loop8.clear();
    loop9.clear();
    loop10.clear();*/
    cout<<"step2\n";
    //cout<<"sum before step 2"<<sum<<"\n";
/*cout<<"JunctionDist before step2:"<<endl;
    for (int i=0;i<NRA;++i){
        for(int j=0;j<NRB;++j){
            cout<<JunctionDist[i][j]<<" ";
        }
        cout<<endl;
    }

    for(size_t i=0;i<NRA;++i){
        cout<<"i"<<i<<endl;
        for (size_t j=0;j<neighA[i].size();++j){
            cout<<neighA[i][j]<<" ";
        }
        cout<<endl;
    }
    */
    //system("pause");
    step=2;

    //cout<<"Conv.size()"<<Conv.size()<<"\n";
    K_forward=1;
    double finalConv2 = KMCconv(conversion,invKeq,K_forward, 0,step);
    cout<<"KMCconv done\n";
    output(step);
    cout<<"Hi\n";
    //FILE * fp;
    //fp=fopen("molecule_size2.csv","a");

    cout<<"c_star\n";
    printf("%.8f\t",c_star);
    cout<<"\n";
    cout<<"loop_frac\n";
    for(size_t l=0;l<loop.size();++l) printf("%.8f\t",loopfrac[l]);
    cout<<"\n";
    cout <<"finalConv2 "<< finalConv2 << "\n";
    //cout << r_thres << "\t";
    cout<< "BoxSize " << BoxSize << "\n";
    //cout << V << "\t";
    //cout << ca0 << "\t";
    cout << "seed " <<seed << "\n";
    cout<< (double)(clock()-c)/CLOCKS_PER_SEC<<endl;
    //_getch();




    //step 3
    cout<<"urA.size()"<<urA.size();
    cout<<"reactedA.size()"<<reactedA.size();

    /*if(urA.size()==0 && reactedA.size()==0){
        cout<<"all sizes are zero";
        system("pause");
    }*/
    Conv.clear();
    time_arr.clear();
    full_reactedA_array.clear();
    all_reactedA_fg.clear();
    Conv_calc.clear();
    Sum.clear();
    
    //loop0.clear();
    loop1.clear();
    loop2.clear();
    loop3.clear();
    loop4.clear();
    /*
    loop5.clear();
    loop6.clear();
    loop7.clear();
    loop8.clear();
    loop9.clear();
    loop10.clear();*/
    cout<<"step3\n";
    //cout<<"sum before step 3"<<sum<<"\n";

    step=3;

    //cout<<"Conv.size()"<<Conv.size()<<"\n";
    K_forward=1;
    double finalConv3 = KMCconv(conversion,0, K_forward,0,step);
    cout<<"KMCconv done\n";
    output(step);
    cout<<"Hi\n";
    //FILE * fp;
    //fp=fopen("molecule_size3.csv","a");

    cout<<"c_star\n";
    printf("%.8f\t",c_star);
    cout<<"\n";
    cout<<"loop_frac\n";
    for(size_t l=0;l<loop.size();++l) printf("%.8f\t",loopfrac[l]);
    cout<<"\n";
    cout <<"finalConv3 "<< finalConv3 << "\n";
    //cout << r_thres << "\t";
    cout<< "BoxSize " << BoxSize << "\n";
    //cout << V << "\t";
    //cout << ca0 << "\t";
    cout << "seed " <<seed << "\n";
    cout<< (double)(clock()-c)/CLOCKS_PER_SEC<<endl;
    //_getch();
    
    //cout<<"total_time"<< (double)(clock()-c)/CLOCKS_PER_SEC<<endl;
    //system("pause");
    return 0;
}

void initialize()
{
    // initialize RNG
    if(DEBUG)
    seed = TEST_SEED;
    else
    seed = time(NULL)+clock()+suffix;
    RNG_initialize(seed);
    // initialize PAB, la
    PAB=(1.0e24/6.02214e23)*pow((3.0/(2.0*3.14159*Nbsq)),1.5);
    la=PAB/ca0;
    //cout<<"la="<<la<<endl;
    //system("pause");
    // initialize box size and rejection related parameters
    V = ((double)NRA*fA / (6.02214e23*ca0)) * 1.0e24; // volume of simulation box in units of nm^3
    BoxSize = pow(V,1.0/3); // Length of simulation box in unit of nm
    double MaxTrials = NRA*fA/10; // not used anywhere else!!
    for(size_t i=0;i<NRA;++i)
    /*
    rA.push_back(vector<double>(3,0.0));// vector with size 3 and all values as 0.0- ie. all NA functional groups have a position given by a vector of length 3
    for(size_t i=0;i<NRA;++i)
        for(size_t x=0;x<3;++x)// 3 coordinates- x,y,and z
            rA[i][x] = rand32()/4294967296.0*BoxSize;// randomly defining the position of each of the crosslinks- where did this formula come from??
    for(size_t i=0;i<NRB;++i) // similarly for B
        rB.push_back(vector<double>(3,0.0));
        for(size_t i=0;i<NRB;++i)
            for(size_t x=0;x<3;++x)
                rB[i][x] = rand32()/4294967296.0*BoxSize;

    */
                // initialize matrix p for determining probability
    p = vector<double>(MAX+1,0);
    if(!LOOPFREE) {
        for(size_t i=1;i<MAX;++i)
            p[i] = 1.0 + la*NRA*fA*pow(i,-1.5);
        p[MAX] = 1.0;
    }         
    p[0] = 1.0;
    // initialize JunctionDist array
    for(size_t i=0;i<NRA;++i){
        JunctionDist.push_back(vector<unsigned short>(NRB,0));
        LoopSize.push_back(vector<unsigned short>(NRB,0));
        
    }
    for(size_t i=0;i<NRB;++i){
        G.push_back(vector<size_t>());
    }
   

   
    // initialize neighA neighB
    for(size_t i=0;i<NRA;++i)
        neighA.push_back(vector<size_t>());
    for(size_t i=0;i<NRB;++i)
        neighB.push_back(vector<size_t>());
        //G[i]=-1;
    // initialize unreacted A and unreacted B vectors
    for(size_t i=0;i<NRA;++i) {
        urA.push_back(i);
    }
    for(size_t i=0;i<NRB;++i) {
        urB.push_back(i);
    }
        // initialize cumulative probability
        sum = NRA*fA * NRB*fB;
        //cout<<"sum initially: "<<sum<<"\n";
        sumA = vector<double>(NRA,fA*NRB*fB);
        
        //cout<<"sum initially: "<<sum<<"\n";
        sumA = vector<double>(NRA,fA*NRB*fB);

        mol = vector<size_t>(NRA+NRB,0);
        largestMol=0;
        for(size_t i=0;i<mol.size();++i) mol[i] = i;
            size_t NNN = NRA*fA>NRB*fB ? NRB*fB : NRA*fA;  // ie. min(NRA*fA,NRB*fB)

        //loop0.reserve(NNN);
        loop1.reserve(NNN);
        loop2.reserve(NNN);
        loop3.reserve(NNN);
        
        loop4.reserve(NNN);
        /*
        loop5.reserve(NNN);
        loop6.reserve(NNN);
        loop7.reserve(NNN);
        loop8.reserve(NNN);
        loop9.reserve(NNN);
        loop10.reserve(NNN);*/
        node1.reserve(NNN);
        node2.reserve(NNN);
        // initialize loop vector
        loop = vector<double>(loopVecSize,0);
        loopfrac = vector<double>(loop.size(),0); // loop.size() should be equal to loopVecSize
        // initialize post KMC calculation variable
}

void output_time_step(size_t ite, size_t step)
{
    FILE *fp;
    char fn[50];
    ////cout<<"writing to file at intermediate time step="<<ite<<"\n";
    if(write_connectivity_data) {
        if(step==1) sprintf(fn,"%s%sKMC_%d.txt",PATH,prefixconnect_1,ite);//,c_star,NRA,NRB,suffix);
        if(step==2) sprintf(fn,"%s%sKMC_%d.txt",PATH,prefixconnect_2,ite);//,c_star,NRA,NRB,suffix);
        if(step==3) sprintf(fn,"%s%sKMC_%d.txt",PATH,prefixconnect_3,ite);//,c_star,NRA,NRB,suffix);
        fp=fopen(fn,"a");
        fprintf(fp,"Node1,Node2\n");
        ////cout<<"node1.size()"<<node1.size()<<"\n";
        for(size_t l=0;l<node1.size();++l){
        //printf("Hello inside loop\n");
        // fprintf(fp,"HEllo");
            fprintf(fp,"%zu,%zu\n",node1[l],node2[l]);
        }
        //printf("Conv[0]",Conv[0],"\n");
        fprintf(fp,"\n");
        fclose(fp);
    }
}
void output(size_t step)
{
    cout<<"wrote to file\n";
    FILE *fp;
    // set file name
    char fn[50]; // no more than 50 characters for file name
    // write loop fraction
    if(writeLoop) {
        sprintf(fn,"%s%scs=%1.4fA%dB%d.csv",PATH,prefixLoopFrac,c_star,NRA,NRB);
        //cout<<fn;
        fp=fopen(fn,"a");
        fprintf(fp,"%.8f,",c_star);
        for(size_t l=0;l<loop.size();++l)
        fprintf(fp,"%.8f,",loopfrac[l]);
        fprintf(fp,"\n");
        fclose(fp);
    }
    // write MW
    if(writeMW) {
        sprintf(fn,"%s%scs=%1.4fA%dB%d_%02d.csv",PATH,prefixMW,c_star,NRA,NRB,suffix);
        fp=fopen(fn,"a");
        fprintf(fp,"Conv,DPw,DPwr,X,SolFrac,BranchFrac,Loop1,Loop2,Loop3,Loop4\n");

        fprintf(fp,"\n");
        fclose(fp);
    }

    if(write_conv_loop) {
        if(step==1)
            sprintf(fn,"%s%scs=%1.4fA%dB%d_%02d.csv",PATH,prefixconv_loop_step1,c_star,NRA,NRB,suffix);
        else if(step==2)
            sprintf(fn,"%s%scs=%1.4fA%dB%d_%02d.csv",PATH,prefixconv_loop_step2,c_star,NRA,NRB,suffix);
        else if(step==3)
            sprintf(fn,"%s%scs=%1.4fA%dB%d_%02d.csv",PATH,prefixconv_loop_step3,c_star,NRA,NRB,suffix);
        fp=fopen(fn,"a");
        //fprintf(fp,"Conv,Loop1,Conv_calc,Sum,loop2,loop3,loop4,loop5,loop6,loop7,loop8,loop9,loop10,loop_sum,loop0\n");
        fprintf(fp,"Conv,Loop1,Conv_calc,Sum,Sum,full_reactedA_array,all_reactedA_fg,loop2, time_arr\n");

        //cout<<"Conv.size()"<<Conv.size()<<"\n";
        //cout<<"Conv_calc.size()"<<Conv_calc.size()<<"\n";
        //cout<<"loop1.size()"<<loop1.size()<<"\n";
        for(size_t l=0;l<loop1.size();++l){

            fprintf(fp,"%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%8f,%8f\n",Conv[l],loop1[l],Conv_calc[l],Sum[l],Sum[l],full_reactedA_array[l],all_reactedA_fg[l],loop2[l],time_arr[l]);
        }
        //printf("Conv[0]",Conv[0],"\n");
        fprintf(fp,"\n");
        fclose(fp);
    }


    if(write_connectivity_data) {

        if(step==1) sprintf(fn,"%s%sKMC.txt",PATH,prefixconnect_1);//,c_star,NRA,NRB,suffix);
        if(step==2) sprintf(fn,"%s%sKMC.txt",PATH,prefixconnect_2);//,c_star,NRA,NRB,suffix);
        if(step==3) sprintf(fn,"%s%sKMC.txt",PATH,prefixconnect_3);//,c_star,NRA,NRB,suffix);
        fp=fopen(fn,"a");
        fprintf(fp,"Node1,Node2\n");
        cout<<"node1.size()"<<node1.size()<<"\n";
        for(size_t l=0;l<node1.size();++l){
        //printf("Hello inside loop\n");
        // fprintf(fp,"HEllo");
            fprintf(fp,"%zu,%zu\n",node1[l],node2[l]);
        }
        //printf("Conv[0]",Conv[0],"\n");
        fprintf(fp,"\n");
        fclose(fp);
    }

}

bool KMCstep_forward()
{
    //cout<<"inside KMCstep_forward"<<endl;
   /* for (size_t A=0;A<NRA;++A){
        for (size_t i=0;i<neighA[A].size();++i){
            size_t B=neighA[A][i];
            if(mol[A]!=mol[B+NRA]){
                cout<<"A"<<A<<"\n";
                cout<<"B"<<neighA[A][i]<<"\n";
                cout<<"mol[A]"<<mol[A]<<"\n";
                cout<<"mol[B]"<<mol[B+NRA]<<"\n";
                cout<<"neighA[A].size()"<<neighA[A].size()<<"\n";
                for (size_t j=0;j<neighA[A].size();++j){
                    cout<<"neighA[A][i]"<<neighA[A][j]<<"\n";}
                cout<<"neighB[B].size()"<<neighB[B].size()<<"\n";
                for (size_t j=0;j<neighB[B].size();++j){
                    cout<<"neighB[B][i]"<<neighB[B][j]<<"\n";}
               // CollectConnected_reverse()
        vector<size_t> AcA_new,BcA_new,AcB_new,BcB_new;
        vector<size_t> dAcA_new,dBcB_new;
    // Collect junctions connected to JA in AcA BcA; connected to JB in AcB BcB
        CollectConnected_reverse(A,B,AcA_new,AcB_new,BcA_new,BcB_new,dAcA_new,dBcB_new); // reverse done- just for testing- to ensure that sum is not updated anywhere
        //cout<<"FOR 144-133:\n";
        vector<size_t> AcA_n,BcA_n,AcB_n,BcB_n;
        vector<size_t> dAcA_n,dBcB_n;
    // Collect junctions connected to JA in AcA BcA; connected to JB in AcB BcB
        CollectConnected_reverse(A,B,AcA_n,AcB_n,BcA_n,BcB_n,dAcA_n,dBcB_n); // reverse done- just for testing- to ensure that sum is not updated anywhere
        
        //cout<<"AcBnew.size()"<<AcB_new.size()<<"\n";
        //for (size_t j=0;j<AcB_new.size();++j){
           // cout<<"AcBnew[i]"<<AcB_new[j]<<"\n";}
        //cout<<"neighB[B].size()"<<neighB[B].size()<<"\n";
        //cout<<"JunctionDist[A][B]"<<JunctionDist[A][B]<<"\n";
        //cout<<"Forward-problem in molecule update somewhere!- case 1";
        mol[A]=mol[B+NRA];
        //cout<<"AFTER CHANGE";
        //cout<<"mol[A]"<<mol[A]<<"\n";
        //cout<<"mol[B]"<<mol[B+NRA]<<"\n";
        //system("pause");
        
        //cout<<"AcBn.size()"<<AcB_n.size()<<"\n";
        //for (size_t j=0;j<AcB_n.size();++j){
            //cout<<"AcBn[i]"<<AcB_n[j]<<"\n";}
        //cout<<"neighB[B].size()"<<neighB[B].size()<<"\n";
        //cout<<"JunctionDist[A][B]"<<JunctionDist[A][B]<<"\n";
        //cout<<"Forward-problem in molecule update somewhere!- case 1";
        
        system("pause");
            }
        }
    }*/

    for (size_t B=0;B<NRB;++B){
        for (size_t i=0;i<neighB[B].size();++i){
            size_t A=neighB[B][i];
            if(mol[A]!=mol[B+NRA]){
                cout<<"B"<<B<<"\n";
                cout<<"A"<<neighB[B][i]<<"\n";
                cout<<"mol[A]"<<mol[A]<<"\n";
                cout<<"mol[B]"<<mol[B+NRA]<<"\n";
                cout<<"neighA[A].size()"<<neighA[A].size()<<"\n";
                cout<<"neighA[A][i]"<<neighA[A][i]<<"\n";
                cout<<"neighB[B].size()"<<neighB[B].size()<<"\n";
                cout<<"JunctionDist[A][B]"<<JunctionDist[A][B]<<"\n";
                
                cout<<"Forward-problem in molecule update somewhere!- case 2\n";
                mol[A]=mol[B+NRA];
                cout<<"AFTER CHANGE";
                cout<<"mol[A]"<<mol[A]<<"\n";
                cout<<"mol[B]"<<mol[B+NRA]<<"\n";
                system("pause");
            }
        }
    }
    // JA JB holds the selected junction for this step
    // idxA idxB holds the index of JA JB in urA urB
    size_t JA,JB,idxA,idxB;
    //size_t numTrials = 0; // not used anyuwhere!!
    // select the junction A and junction B to react
    SelectJunct_forward(JA,JB,idxA,idxB); // JA,JB,idxA,idxB are passed as pointers- so, these are not initialized earlier
    UpdateConnectivityData_forward(JA,JB);
    //cout<<"Forward"<<"JA"<<JA<<"JB"<<JB<<"\n";
    // updating adjacency list
    // this is relevant only for bifunctional polymer chains. Will have to modify slightly for other functionalities
    if(neighA[JA].size()>0){
        size_t JB1=neighA[JA][0]; // since this is before nieghbour update, only the first neighbour has to be considered
        // only valid for bifunctional polymers 
        //cout<<"calculated JB1="<<JB1<<endl;
        //cout<<"JB:"<<JB<<endl;
        //cout<<"G[0].size()"<<G[0].size()<<endl;
       // cout<<"G[1].size()"<<G[1].size()<<endl;
        //cout<<"G[2].size()"<<G[2].size()<<endl;

        //cout<<"G[JB].size()"<<G[JB].size()<<endl;
        G[JB].push_back(JB1);
       // cout<<"updated G[JB]"<<endl;
        //cout<<"G[JB1].size()"<<G[JB1].size()<<endl;
        G[JB1].push_back(JB);
        //cout<<"G[JB1].size()"<<G[JB1].size()<<endl;
        //cout<<"updated G[JB1]"<<endl;
    }
    //cout<<"after G update in KMCstep_forward"<<endl;
    
   
    
      


    // Updates the sum of relative probabilities of unreacted A-B pairs,
    //by deleting the probabilities for the A and B reacting this step
    //cout<<"JA"<<JA<<"\n";
    //cout<<"JB"<<JB<<"\n";
    UpdateSum_forward(JA,JB,idxA,idxB);
    //cout<<urA.size()<<"\n";
    // Update Loop information
    UpdateLoop_forward(JA,JB);
    // Find all junctions connected to JA and JB
    // AcA holds all A junctions that is connected to JA
    // BcA holds all B junctions that is connected to JA
    // AcB holds all A junctions that is connected to JB
    // BcB holds all B junctions that is connected to JB
    vector<size_t> AcA,BcA,AcB,BcB;
    vector<size_t> dAcA,dBcB;
    // Collect junctions connected to JA in AcA BcA; connected to JB in AcB BcB
    CollectConnected(JA,JB,AcA,AcB,BcA,BcB,dAcA,dBcB);
    // Updates connectivity by updating neighbor list of JA and JB
    neighA[JA].push_back(JB);
    neighB[JB].push_back(JA);
    bond new_bond;
    new_bond.JA=JA;
    new_bond.JB=JB;
    all_bonds.push_back(new_bond);
    //cout<<"update neighbors for JA and JB: "<<JA<<"  "<<JB<<"\n";
    //cout<<"neighA[JA].size()"<<neighA[JA].size()<<"\n";
    //cout<<"neighB[JB].size()"<<neighB[JB].size()<<"\n";
    // update JunctionDist
    UpdateJuncDist_forward(JA,JB,BcA,AcB); // same function definition, but different inputs- need to decide which one is 
    // being used based on the inputs only
    UpdateJuncDist_common(JA,JB,AcA,BcB,dAcA,dBcB);
    // Update the molecule grouping info and molecule sizes
   // cout<<"before molecule update\n";
    //cout<<"mol[JA]"<<mol[JA]<<"\n";
    //cout<<"mol[JB]"<<mol[JB+NRA]<<"\n";
    //size_t mA1=mol[JA];
    //size_t mB1=mol[JB+NRA];
    UpdateMol_forward(JA,JB,AcA,AcB,BcA,BcB);
    //cout<<"after molecule update\n";
    //cout<<"mol[JA]"<<mol[JA]<<"\n";
    //cout<<"mol[JB]"<<mol[JB+NRA]<<"\n";
    //size_t mA2=mol[JA];
    //size_t mB2=mol[JB+NRA];
   /* if((mA1!=mB1) && (mA1==mA2) && (mB1==mB2)){
        cout<<"PROBLEM in moleucle update\n";
        system("pause");
    }*/
    

    /*cout<<"AcA[j]:";
        for(size_t j=0;j<AcA.size();++j) cout<<AcA[j]<<" ";
        cout<<endl;
    cout<<"BcA[j]:";
        for(size_t j=0;j<BcA.size();++j) cout<<BcA[j]<<" ";
        cout<<endl;
    cout<<"AcB[j]:";
    for(size_t j=0;j<AcB.size();++j) cout<<AcB[j]<<" ";
        cout<<endl;
    cout<<"BcB[j]:";
        for(size_t j=0;j<BcB.size();++j) cout<<BcB[j]<<" ";
        cout<<endl;
    cout<<"JunctionDist[JA][JB]"<<JunctionDist[JA][JB]<<endl;
    cout<<"mol[101]"<<mol[101]<<endl;
    cout<<"mol[43]"<<mol[43]<<endl;
        //system("pause");
    */
     /*for(size_t i=0;i<reactedA.size();++i){
        if(neighA[reactedA[i]].size()<fA){
            cout<<"in forward neighA[reactedA[i]].size<fA\n";
            cout<<"in forward reactedA[i]"<<reactedA[i]<<"\n";
            cout<<"neighA[reactedA[i]].size()"<<neighA[reactedA[i]].size()<<endl;
            system("pause");
        }
    }*/

    /*cout<<"JunctionDist_forward:"<<endl;
    for (int i=0;i<NRA;++i){
        for(int j=0;j<NRB;++j){
            cout<<JunctionDist[i][j]<<" ";
        }
        cout<<endl;
    }*/
   //cout<<"sum after forward"<<sum<<endl;
    return true;
}






bool KMCstep_reverse()
{
    //cout<<"inside KMCstep_reverse"<<endl;
    //sumrA_check--;
    // select junctions
    // JA JB holds the selected junction for this step
    // idxA idxB holds the index of JA JB in urA urB

    if(all_bonds.size()>0){  // reverse step can run only if there is at least one bond in the system, otherwise- cannot run
    size_t JA,JB,idxA,idxB;
    JA=MAX;
    JB=MAX;
    for (size_t A=0;A<NRA;++A){
        for (size_t i=0;i<neighA[A].size();++i){
            size_t B=neighA[A][i];
            if(mol[A]!=mol[B+NRA]){
                cout<<"A"<<A<<"\n";
                cout<<"B"<<neighA[A][i]<<"\n";
                cout<<"mol[A]"<<mol[A]<<"\n";
                cout<<"mol[B]"<<mol[B+NRA]<<"\n";
                cout<<"neighA[A].size()"<<neighA[A].size()<<"\n";
                cout<<"neighA[A]:"<<endl;
                for(size_t i=0;i<neighA[A].size();++i){
                    cout<<neighA[A][i]<<" ";
                }
                cout<<endl;
                cout<<"neighB[B]:"<<endl;
                for(size_t i=0;i<neighB[B].size();++i){
                    cout<<neighB[B][i]<<" ";
                }
                cout<<endl;
                cout<<"neighA[A][0]"<<neighA[A][0]<<"\n";
                cout<<"neighB[B].size()"<<neighB[B].size()<<"\n";
                cout<<"JunctionDist[A][B]"<<JunctionDist[A][B]<<"\n";
                cout<<"Reverse-problem in molecule update somewhere!-case 1\n";
                mol[A]=mol[B+NRA];
                cout<<"AFTER CHANGE";
                cout<<"mol[A]"<<mol[A]<<"\n";
                cout<<"mol[B]"<<mol[B+NRA]<<"\n";
                system("pause");
            }
        }
    }

        
    /*for (size_t B=0;B<NRB;++B){
        for (size_t i=0;i<neighB[B].size();++i){
            size_t A=neighB[B][i];
            if(mol[A]!=mol[B+NRA]){
                cout<<"B"<<B<<"\n";
                cout<<"A"<<neighB[B][i]<<"\n";
                cout<<"mol[A]"<<mol[A]<<"\n";
                cout<<"mol[B]"<<mol[B+NRA]<<"\n";
                cout<<"neighA[A].size()"<<neighA[A].size()<<"\n";
                cout<<"neighA[A][i]"<<neighA[A][i]<<"\n";
                cout<<"neighB[B].size()"<<neighB[B].size()<<"\n";
                cout<<"JunctionDist[A][B]"<<JunctionDist[A][B]<<"\n";
                
                cout<<"Reverse-problem in molecule update somewhere!- case 2\n";
                mol[A]=mol[B+NRA];
                cout<<"AFTER CHANGE";
                cout<<"mol[A]"<<mol[A]<<"\n";
                cout<<"mol[B]"<<mol[B+NRA]<<"\n";
                system("pause");
            }
        }
    }*/

    /*for(size_t i=0;i<reactedA.size();++i){
        if(neighA[reactedA[i]].size()<fA){
            cout<<"neighA[reactedA[i]].size<fA\n";
            cout<<"reactedA[i]"<<reactedA[i]<<"\n";
            system("pause");
        }
    }

    for(size_t i=0;i<reactedB.size();++i){
        if(neighB[reactedB[i]].size()==0){
            cout<<"neighB[reactedB[i]].size=0\n";
            cout<<"i"<<i<<"\n";
            system("pause");
        }
    }*/

    // select the junction A and junction B to react
    SelectJunct_reverse(JA,JB);//,idxA,idxB); // JA,JB,idxA,idxB are passed as pointers- so, these are not initialized earlier

    ////cout<<"Reverse:"<<"JA"<<JA<<"JB"<<JB<<"\n";

    //IMPORTANT!!:  NEED TO PUT THIS FOLLOWING JUNCTION SELECTION STEP IN SelectJunct_reverse
    //cout<<"JA"<<JA<<"\n";
    //cout<<"JB"<<JB<<"\n";

    //G.erase(std::remove(reactedA.begin(),reactedA.end(), JA),reactedA.end());
    // doing this before updating neighbors

     // updating adjacency list
    // this is relevant only for bifunctional polymer chains. Will have to modify slightly for other functionalities
    if(neighA[JA].size()==fA){
        size_t JB1=neighA[JA][0];
        size_t JB2=neighA[JA][1];
        size_t JB_bar=(JB2==JB) ? JB1 : JB2; // the one which is not JB

        for (size_t i=0;i<G[JB].size();++i){
                if(G[JB][i]==JB_bar){
      
                    G[JB].erase(G[JB].begin()+i);
                    break;
                }
            }
            for (size_t i=0;i<G[JB_bar].size();++i){
                if(G[JB_bar][i]==JB){
          
                    G[JB_bar].erase(G[JB_bar].begin()+i);
                    break;
                }
            }

/*        
        if(JB2==JB){

            for (size_t i=0;i<G[JB].size();++i){
                if(G[JB][i]==JB1){
            //cout<<"reverse JA"<<JA<<"\n";
                    G[JB].erase(G[JB].begin()+i);
                    break;
                }
            }
            for (size_t i=0;i<G[JB1].size();++i){
                if(G[JB1][i]==JB){
            //cout<<"reverse JA"<<JA<<"\n";
                    G[JB1].erase(G[JB1].begin()+i);
                    break;
                }
            }

           // G[JB].erase(std::remove(G[JB].begin(),G[JB].end(), JB1),G[JB].end());
            //cout<<"erased JB1 from JB"<<endl;
            //G[JB1].erase(std::remove(G[JB1].begin(),G[JB1].end(), JB),G[JB1].end());
            //cout<<"erased JB from JB1"<<endl;
           // G[JB].remove(JB1);
            //G[JB1].remove(JB);
        }
        else{

            // DONT USE ERASE REMOVE IDIOM HERE- WILL DELETE ALL INSTANCES OF THE ELEMENT!!!

            for (size_t i=0;i<G[JB].size();++i){
                if(G[JB][i]==JB2){
            //cout<<"reverse JA"<<JA<<"\n";
                    G[JB].erase(G[JB].begin()+i);
                    break;
                }
            }
            for (size_t i=0;i<G[JB2].size();++i){
                if(G[JB2][i]==JB){
            //cout<<"reverse JA"<<JA<<"\n";
                    G[JB2].erase(G[JB2].begin()+i);
                    break;
                }
            }
            //G[JB].erase(std::remove(G[JB].begin(),G[JB].end(), JB2),G[JB].end());
            //cout<<"erased JB2 from JB"<<endl;
            //G[JB2].erase(std::remove(G[JB2].begin(),G[JB2].end(), JB),G[JB2].end());
            //cout<<"erased JB from JB2"<<endl;
            //G[JB].remove(JB2);
            //G[JB2].remove(JB);
        }
*/


        //cout<<"after G update in KMCstep_reverse"<<endl;
        
    }

    if(JA==MAX || JB==MAX){
        cout<<"JA: "<<JA;
        cout<<"JB: "<< JB;
        cout<<"either JA or JB or both are MAX\n";
        system("pause");
        return false;
    }
    else{
        //cout<<"checking molecule before update junction reverse, will pause after this pause if found wrong!\n";
        /*if(mol[JA]!=mol[JB+NRA]){
            cout<<"problem in molecule similarity before update junction reverse\n";
            system("pause");
        }*/


    // Updates the sum of relative probabilities of unreacted A-B pairs,
    //by deleting the probabilities for the A and B reacting this step
    //size_t idx_check = count(urA.begin(), urA.end(), JA); 
    //if(idx_check<=0) urA.push_back(JA); // implies not found 
    
    if(neighA[JA].size() == fA) {
        urA.push_back(JA); // this is equivalent. because if the functionality is less than fA, then JA would already havebeen part of urA
        reactedA.erase(std::remove(reactedA.begin(),reactedA.end(), JA),reactedA.end());
        
        //urA.erase(urA.begin()+idxA);
    }

    //idx_check = count(urB.begin(), urB.end(), JB); 
    //if(idx_check<=0) urB.push_back(JB); // implies not found
    if(neighB[JB].size() == fB) {
        urB.push_back(JB);
        reactedB.erase(std::remove(reactedB.begin(),reactedB.end(), JB),reactedB.end());
        //urB.push_back(JB);
        //reactedB.push_back(JB);
        //urB.erase(urB.begin()+idxB);
    }
    sort(urA.begin(), urA.end());
    sort(urB.begin(), urB.end());


    
    // Update Loop information
    UpdateLoop_reverse(JA,JB);
    // Find all junctions connected to JA and JB
    // AcA holds all A junctions that is connected to JA
    // BcA holds all B junctions that is connected to JA
    // AcB holds all A junctions that is connected to JB
    // BcB holds all B junctions that is connected to JB
    vector<size_t> AcA,BcA,AcB,BcB;
    vector<size_t> dAcA,dBcB;

    for (size_t i=0;i<neighA[JA].size();++i){
        if(neighA[JA][i]==JB){
            //cout<<"reverse JB"<<JB<<"\n";
            neighA[JA].erase(neighA[JA].begin()+i);
            //cout<<"erased JB from JA connection: JB="<<JB<<" JA="<<JA<<endl;
            //cout<<"neighA[JA].size()="<<neighA[JA].size()<<endl;
            //if(JA==185){cout<<"neighA[185] deleted in reverse"<<endl; system("pause");}
            break;// this break is required to make sure that only one bond is erased!!- EG. IN CASE OF LOOP
        }
    }
    bool loop_break=0;
    for (size_t i=0;i<neighA[JA].size();++i){
        if(neighA[JA][i]==JB) {//JunctionDist[JA][JB]=1;// since they are neighbours, junctions are NOT being temporariliy disconnected!!
        //cout<<"JA and JB are neighbours, hence not temporarily disconnected!!\n";
        //cout<<"LOOP- neighbors even after disconnected\n";
        loop_break=1;

        break;
        }
    }

    for (size_t i=0;i<neighB[JB].size();++i){
        if(neighB[JB][i]==JA){
            //cout<<"reverse JA"<<JA<<"\n";
            neighB[JB].erase(neighB[JB].begin()+i);
            break;
        }
    }

    /*for(size_t i=0;i<reactedA.size();++i){
        if(neighA[reactedA[i]].size()<fA){
           // cout<<"after reverse neighA[reactedA[i]].size<fA\n";
           // cout<<"after reverse reactedA[i]"<<reactedA[i]<<"\n";
           // cout<<"reactedA[i]"<<reactedA[i]<<endl;
            system("pause");
        }
    }*/




    // Collect junctions connected to JA in AcA BcA; connected to JB in AcB BcB
    //vector<size_t> AcA_rev,BcA_rev,AcB_rev,BcB_rev;
    //vector<size_t> dAcA_rev,dBcB_rev;

    //cout<<"Before CollectConnected_reverse-1\n";
    //CollectConnected_reverse(JA,JB,AcA_rev,AcB_rev,BcA_rev,BcB_rev,dAcA_rev,dBcB_rev);
    //cout<<"Done CollectConnected_reverse-1\n";

    // Updates connectivity by updating neighbor list of JA and JB
    // after reverse rxn- JA and JB are no longer connected!!

    //IMPORTANT: DONT USE ERASE REMOVEIDOM HERE- BECAUSE IT IS ERASING ALL INSTANCES OF THE SPECIFIC JUCNTION- BUT
    //PHYSICALLY THAT IS NOT CORRECT
    // HAVING MULTIPLE INSTANCES MEANS THAT THERE ARE MORE THAN ONE CONNECTIONS- AND ERASING ALL OF THEM WOULD MEAN THAT 
    // I AM DISCONNECTING ALL CONNECTIONS- WHICH IS WRONG!!!
    //neighA[JA].erase(std::remove(neighA[JA].begin(),neighA[JA].end(), JB),neighA[JA].end());
    //neighB[JB].erase(std::remove(neighB[JB].begin(),neighB[JB].end(), JA),neighB[JB].end());

 
    //cout<<"Before UpdateJuncDist_reverse-1\n";
    UpdateConnectivityData_reverse(JA,JB);
    UpdateJuncDist_reverse(JA,JB);//,BcA_rev);//,AcB_rev,neighA,neighB); // same function definition, but different inputs- need to decide which one is 
    // being used based on the inputs only
    //cout<<"After UpdateJuncDist_reverse-1\n";
  

    // Update the molecule grouping info and molecule sizes
    //cout<<"JA:"<<JA<<"\n";
    //cout<<"JB:"<<JB<<"\n";
    //cout<<"JunctionDist[0][0]"<<JunctionDist[0][0]<<endl;
    vector<size_t> AcA_new,BcA_new,AcB_new,BcB_new;
    vector<size_t> dAcA_new,dBcB_new;
    // Collect junctions connected to JA in AcA BcA; connected to JB in AcB BcB
    CollectConnected_reverse(JA,JB,AcA_new,AcB_new,BcA_new,BcB_new,dAcA_new,dBcB_new);
    // Doing this again so that collect connected is implemented after junctiondistance update!!
    // Update the molecule grouping info and molecule sizes

    UpdateMol_reverse(JA,JB,AcA_new,AcB_new,BcA_new,BcB_new);
    UpdateSum_reverse(JA,JB);
    //cout<<"sum after reverse"<<sum<<endl;
 
    return true;
    }

    }
    
    else{
        return false;
    }

}


double KMCconv(double conversion, double invKeq, double K_forward, double K, size_t step)
{
    double time_curr=0;
    // do KMCstep until designated conversion
    //NA=NRA*fA,NB=NRB*fB;
    double End = (1-conversion)*start;
    if(LOOPFREE) End = 1 + start - ( start/fA+start/fB );
    
    double currConv=1.0;
    double currConv_calc=1.0;
    size_t totalbond = start - ceil(End);
    //size_t i_rxn=start;// index for tracking conversion
    size_t i_temp;
    
    size_t i_start=start;
    size_t count_ite=0;
    if(invKeq>0) i_start=2*start;
    for( size_t i = i_start; i > End; i--) {

    /*
    double sumfA=0;//sum of all reacted B "functional groups"- useful for forward rxns
    double sumrA=0;// sum of all reacted A "functional groups"- useful for reverse rxns
    for(size_t k=0; k<reactedA.size(); ++k){
        size_t num_reactedA=neighA[reactedA[k]].size();
        sumrA=sumrA+num_reactedA;
    }
    for(size_t k=0; k<urA.size(); ++k){
        size_t num_reactedA=neighA[urA[k]].size();
        size_t num_unreactedA=fA-num_reactedA;
        sumfA=sumfA+num_unreactedA;
        sumrA=sumrA+num_reactedA;
    }
    double sumfB=0;
    double sumrB=0;// sum of all reacted B "functional groups"
    for(size_t k=0; k<reactedB.size(); ++k){
        size_t num_reactedB=neighB[reactedB[k]].size();
        sumrB=sumrB+num_reactedB;
    }
    for(size_t k=0; k<urB.size(); ++k){
        
        size_t num_reactedB=neighB[urB[k]].size();
        size_t num_unreactedB=fB-num_reactedB;
        sumfB=sumfB+num_unreactedB;
        sumrB=sumrB+num_reactedB;
    }

    if(sumrA!=sumrB) {
        cout<<"PROBLEM!! sumrA!=sumrB"<<"\n";
        cout<<"sumrA: "<<sumrA<<"\n";
        cout<<"sumrB "<<sumrB<<"\n";
        cout<<"sumfA: "<<sumfA<<"\n";
        cout<<"sumfB "<<sumfB<<"\n";
        system("pause");
    }
    */

    double sum_f=sum;///((6.02214e23*V*1.0e-24)*(6.02214e23*V*1.0e-24));/// concentration is in moles/m3//- THIS SHOULD BE EQUAL TO THE SUM (Sum of probabilities of unreacted junctions)- AND NOT MERELY THE NUMBER OF BONDS AVAILABLE FOR REACTION
    
    double sum_r=all_bonds.size()*(6.02214e23*V*1.0e-24);//sumrA; // // THIS SHOULD BE EQUAL TO THE NUMBER OF AVAILABLE BONDS!! - NOT THIS!!

    double sum_total=sum_f+sum_r*invKeq;
    double rnd_num_tmp=rand32()/4294967296.0;
    time_curr=time_curr+log(1/rnd_num_tmp)*(NRA*fA*NRB*fB/sum_total);

    //cout<<"time_curr="<<time_curr<<endl;
    //if(sum_f==0.0 && invKeq==0.0) {system("pause"); break; }
    double P_f;
    double P_a;
    //if(invKeq==0.0) {P_f=1.0; cout<<"P_f=1 in first loop\n";}
    //else{
   /* cout<<"sum"<<sum<<endl;

    double sum_of_elems = std::accumulate(vector.begin(), vector.end(), 0.0);
    cout<<"sum(sumA)"<<sum_of_elems<<endl;*/
    if(invKeq==0.0) { // only dissoc forward
        P_f=((double) sum_f/(sum_f+0)<1.0? (double) sum_f/(sum_f+0) : 1.0);
        //cout<<"Step 1: sum_f/sum_r"<<sum_f/(sum_r)<<" conversion_calc="<<currConv_calc<<" conversion="<<currConv<<endl;
        //system("pause");
        //min((double) sum_f/(sum_f+sum_a), 1.0);
       // cout<<"Pf"<<P_f;
       // system("pause");
        P_a=1-P_f;
    }
    else if(sum_f==0.0) { // 
        P_f=0.0;// sum_f=0 means that no unreacted left, all reacted, so we must do reverse rxn only!!
        P_a=min((double)0/(0+invKeq*sum_r),1.0);
    // this is just for the first reverse rxn step!!
    }
    else {
    // P_f=min((double) 1.0/(1.0+invKeq*sum_r*V/sum_f), 1.0);//min(1.0/(1.0+invKeq*sum_r*V/sum_f),1.0)
        P_f=K_forward*min((double) sum_f/(sum_f+0+ invKeq*sum_r), 1.0);
        //cout<<"Step 2: sum_f/sum_r*invKeq"<<sum_f/(sum_r*invKeq)<<" conversion="<<currConv_calc<<" conversion="<<currConv<<endl;
        //system("pause");
        P_a=min((double) 0/(sum_f+0+ invKeq*sum_r), 1.0);
        //cout<<"P_f"<<P_f<<endl;
        //if((P_f<0.51) && (P_f>0.49)){ // this means that rate of both forward and reverse reactions is close to equal
        
        
    }
    //P_a=0;
    double rand_num=rand32()/4294967296.0;
    //cout<<"P_f"<<P_f<<"\n";
    //cout<<"P_a"<<P_a<<"\n";
    //cout<<"rand_num"<<rand_num<<"\n";

    // sum calculation using naive method:
    
    /*
    double sum_test_common=0;

    vector<double> sumA_test_common; 
    sumA_test_common = vector<double>(NRA,0);
    for(size_t j=0;j<urA.size();++j){
        size_t JA_loop=urA[j];
        sumA_test_common[JA_loop]=0;
        double probAB;
        for(size_t j=0;j<urB.size();++j){
            size_t JB_loop=urB[j];
            probAB=p[dist(JunctionDist[JA_loop][JB_loop])]*(fA-(double)neighA[JA_loop].size())*(fB-(double)neighB[JB_loop].size());
            sum_test_common+=probAB;
            sumA_test_common[JA_loop]+=probAB;
        }
    }
    for(size_t j=0;j<reactedA.size();++j){
        sumA_test_common[reactedA[j]]=0;
    }*/

    if(rand_num<P_f){
        // forward rxn- dissoc

    //cout<<"rxn_1\n";
    //cout<<"sum forward "<<sum<<"\n";

    if(KMCstep_forward()) {
        //cout<<"KMCstep_forward() done"<<endl;
        //cout<<"i_rxn"<<i_rxn<<endl;
        i_rxn--;
        currConv = (double)(start-i_rxn+1)/(double)start;
            

    //currConv_calc= 1.0-(double) urA.size()/(NRA);
    double num_ur_A_junc=0;
    for (size_t j=0;j<urA.size();++j){
            num_ur_A_junc+=fA-neighA[urA[j]].size();
    }
    currConv_calc= 1.0-(double) num_ur_A_junc/((double)(NA));
   
    //cout<<"i_rxn: "<<i_rxn<<"\n";
    //cout<<"currConv"<<currConv<<"\n";
    //cout<<"currConv_calc"<<currConv_calc<<"\n";
    //cout<<"forward: num_ur_A_junc "<<num_ur_A_junc;
    //cout<<" (start-i_rxn+1)"<<(i_rxn)<<" urA.size()"<<urA.size()<<endl;
    if(!((start-i+1)%MWfreq)) {//cout<<"writing to file at this i: Rxn1: "<<i<<"\n";
        updateWriteData(currConv,currConv_calc, totalbond, time_curr, invKeq);// basically means that do this step when (start-i+1)%MWfreq) is 0
    }
    if(sum < 1 ) {cout<<"breaking KMCConv -forward because sum<1\n";cout<<"sum:"<<sum<<"\n";break;}
    } 
    else { cout<<"KMCstep_forward is false"; break; }
    }

    
    else{ // reverse rxn
        //cout<<"reverse rxn\n";
    

    //doing reverse reaction
    if(KMCstep_reverse()) {
        //cout<<"KMCstep_reverse done"<<endl;
        i_rxn++;
        currConv = (double)(start-i_rxn+1)/(double)start;
            
        //currConv_calc= 1.0-(double) urA.size()/(NRA);
        double num_ur_A_junc=0;
        for (size_t j=0;j<urA.size();++j){
            num_ur_A_junc+=fA-neighA[urA[j]].size();
        }

        /*
        double num_ur_A_junc=0;
        for (size_t j=0;j<urA.size();++j){
            cout<<"urA[j]"<<urA[j]<<endl;
            cout<<"fA="<<fA<<endl;
            cout<<"neighA[urA[j]].size()"<<neighA[urA[j]].size()<<endl;
            for(size_t k=0; k<neighA[urA[j]].size();++k){
                cout<<"k="<< k ;
                cout<<" neighA[urA[j]][k]"<<neighA[urA[j]][k]<<endl;
            }
            if(neighA[urA[j]].size()==0){// means unreacted - NO NEIGHBORS
                //num_ur_A_junc+=fA;
                num_ur_A_junc+=fA-neighA[urA[j]].size();
            }
        }*/
        currConv_calc= 1.0-(double) num_ur_A_junc/((double)(NRA*fA));
        //cout<<"reverse: num_ur_A_junc "<<num_ur_A_junc;
        //cout<<" (start-i_rxn+1)"<<(i_rxn)<<" urA.size()"<<urA.size()<<endl;
        //system("pause");
        //cout<<"i_rxn: "<<i_rxn<<"\n";
        //cout<<"currConv"<<currConv<<"\n";
        //cout<<"currConv_calc"<<currConv_calc<<"\n";
        //if(!((start-i+1)%MWfreq)) updateWriteData(currConv,currConv_calc, totalbond);// basically means that do this step when (start-i+1)%MWfreq) is 0
        if(!((start-i+1)%MWfreq)) {//cout<<"writing to file at this i: Rxn2: "<<i<<"\n";
            updateWriteData(currConv,currConv_calc, totalbond, time_curr, invKeq);// basically means that do this step when (start-i+1)%MWfreq) is 0
        }
        // if this consition is also there for forward- then this will never be encountered in reverse!
    } 
    else { cout<<"KMCstep_reverse is false, mostly because all bonds have been broken\n"; break; }
    }
    //cout<<endl;
     /*if(step==1 && currConv_calc>0.99){cout<<endl;cout<<"step: "<<step<<endl;cout<<"(i_start-i): "<<(i_start-i)<<endl;cout<<"i_start: "<<i_start<<endl;//system("pause");
     break;}
     if(step==2 && currConv_calc<0.01){cout<<endl;cout<<"step: "<<step<<endl;cout<<"(i_start-i): "<<(i_start-i)<<endl;cout<<"i_start: "<<i_start<<endl;//system("pause");
     break;} // have to chose this conversion correctly
      if(step==3 && currConv_calc>0.99){cout<<endl;cout<<"step: "<<step<<endl;cout<<"(i_start-i): "<<(i_start-i)<<endl;cout<<"i_start: "<<i_start<<endl;//system("pause");
      break;}*/
    i_temp=i;
    //cout<<"currConv_calc"<<currConv_calc<<"\n";
    
    if(count_ite % wrt_step==0){
        output_time_step(count_ite,step); // write network to file
    }
    count_ite=count_ite+1;
    }
    //system("pause");
    cout<<"total number of rxn steps: i_start-i_temp "<<i_start-i_temp<<endl;
    cout<<"step: "<<step<<endl;cout<<"i_temp: "<<i_temp<<endl;cout<<"i_start: "<<i_start<<endl;
    cout<<"currConv_calc"<<currConv_calc<<endl;
    // calculate loop fraction from loop[]
    cout<<"all rxn steps done\n";
    for(size_t l=0;l<loop.size();++l)
        loopfrac[l]=((double)loop[l])/totalbond;
    return currConv;
}

void UpdateConnectivityData_forward(size_t &JA,size_t &JB) //TODO
{
char fn[50];
//sprintf(fn,"network.csv");//,PATH,"NW_",c_star,currConv,NRA,NRB,suffix);
//FILE *fp;
//fp=fopen(fn,"a");

if(neighA[JA].size()>0){// means there is at least one B4 connected to A2
    size_t n1=neighA[JA][0] ; // node1-there can be at most one B node connected to A
    size_t n2=JB; // node2- to be connected
    //fprintf(fp,"%d,%d\n",n1,n2);
    node1.push_back(n1);
    node2.push_back(n2);
    //cout<<"updated to list "<<n1<<" "<<n2<<"\n";
    //system("pause");
}
// if neighA[JA].size()==0 then this means that a free A node is just getting connected with B, and hence, no change in connectivity pattern of B
if(neighA[JA].size()==2){

    cout<<"PROBLEM in junction selection\n";
    system("pause");
}
/*
for(size_t i=0;i<neighA.size();++i) {
fprintf(fp,"A%d",i);
for(size_t j=0;j<neighA[i].size();++j) {
fprintf(fp,",B%d",neighA[i][j]);
}
}*/
}

void UpdateConnectivityData_reverse(size_t &JA,size_t &JB) //TODO
{

// this should be done after updating neighA and neighB
char fn[50];
//sprintf(fn,"network.csv");//,PATH,"NW_",c_star,currConv,NRA,NRB,suffix);
//FILE *fp;
//fp=fopen(fn,"a");

if(neighA[JA].size()>0){// means there is at least one B4 connected to A2
    size_t n1=neighA[JA][0] ; // node1-there can be at most one B node connected to A
    // this is the other node connected to JA- which is still attached
    // this node should now get disconnected from JB
    size_t n2=JB; // node2- to be connected
    //fprintf(fp,"%d,%d\n",n1,n2);

    for (size_t t = 0; t < node1.size(); ++t) {
        if ((node1[t] == n1 && node2[t] == n2) || (node1[t] == n2 && node2[t] == n1)) {
            node1.erase(node1.begin() + t);
            node2.erase(node2.begin() + t);
            break;  // remove only the first matching pair
        }
    }
    //cout<<"updated to list "<<n1<<" "<<n2<<"\n";
    //system("pause");
}
// if neighA[JA].size()==0 then this means that a free A node is just getting connected with B, and hence, no change in connectivity pattern of B
if(neighA[JA].size()==2){

    cout<<"PROBLEM in junction selection\n";
    system("pause");
}
/*
for(size_t i=0;i<neighA.size();++i) {
fprintf(fp,"A%d",i);
for(size_t j=0;j<neighA[i].size();++j) {
fprintf(fp,",B%d",neighA[i][j]);
}
}*/
}


void updateWriteData(double currConv,double currConv_calc, size_t totalbond, double time_curr, double invKeq)
{
    //size_t NA=NRA*fA,NB=NRB*fB;
    double totalWW=0,totalW=0;

    //loop0.push_back(loop[0]/totalbond);
    loop1.push_back(loop[1]/totalbond);
    loop2.push_back(loop[2]/totalbond);
    loop3.push_back(loop[3]/totalbond);
    loop4.push_back(loop[4]/totalbond);
    
    /*
    loop5.push_back(loop[5]/totalbond);
    loop6.push_back(loop[6]/totalbond);
    loop7.push_back(loop[7]/totalbond);
    loop8.push_back(loop[8]/totalbond);
    loop9.push_back(loop[9]/totalbond);
    loop10.push_back(loop[10]/totalbond);*/

    //cout<<"loop1 value"<<loop[1]/totalbond<<"\n";
    Conv.push_back(currConv);

    time_arr.push_back(time_curr);
    full_reactedA_array.push_back(reactedA.size());
    Conv_calc.push_back(currConv_calc);
    double sum_f_=sum;///((6.02214e23*V*1.0e-24)*(6.02214e23*V*1.0e-24));/// concentration is in moles/m3//- THIS SHOULD BE EQUAL TO THE SUM (Sum of probabilities of unreacted junctions)- AND NOT MERELY THE NUMBER OF BONDS AVAILABLE FOR REACTION
    
    double sum_r_=all_bonds.size()*(6.02214e23*V*1.0e-24);//sumrA; // // THIS SHOULD BE EQUAL TO THE NUMBER OF AVAILABLE BONDS!! - NOT THIS!!
   // cout<<"sum_f_/sum_r_"<<sum_f_/sum_r_<<endl;

    Sum.push_back(sum_f_+sum_r_*invKeq);
   

    double num_react_A=0; // number of reacted groups
    for (size_t j=0;j<urA.size();++j){
            num_react_A+=neighA[urA[j]].size();
    }
    num_react_A+=fA*reactedA.size();
    all_reactedA_fg.push_back(num_react_A);
}

void SelectJunct_forward(size_t &JA,size_t &JB,size_t &idxA,size_t &idxB)
{
    double stop = rand32()/4294967296.0*sum; //rand32()/4294967296.0 is a uniformly distributed random number between 0 and 1
    //stop = 0.5*sum;// this is only for debugging- delete this line later!!
    double cump = 0;
    int JA_found=0;
    int JB_found=0;
    for(size_t i=0;i<urA.size();++i) {
        if(cump+sumA[urA[i]] >= stop) {
        // cout<<"JA done: sumA[urA[i]]="<<sumA[urA[i]]<<"\n";
        JA = urA[i];
        idxA = i;
        JA_found=1;
        //cout<<"JA= "<<JA<<endl;
        break;
        } else {
        cump += sumA[urA[i]];
        //cout<<"sumA[urA[i]]"<<sumA[urA[i]]<<"\n";

        }
    }
    if(JA_found==0){
        cout<<"JA_found=0\n";
        cout<<"stop"<<stop<<"\n";
        cout<<"cump"<<cump<<"\n";
        cout<<"urA.size()"<<urA.size()<<"\n";
        system("pause");
    }
    for(size_t i=0;i<urB.size();++i) {
        JB = urB[i];
        double tmpP = p[dist(JunctionDist[JA][JB])] * (double)(fB - neighB[JB].size()) * (double)(fA -neighA[JA].size());// this is degeneracy
        if(cump+tmpP >= stop) {
        idxB = i;
        JB_found=1;

        //cout<<"i "<<i<<endl;
        //cout<<"tmpP = "<<tmpP<<"\n";


        break;
        } else {
        cump += tmpP;
        
        }
        
    }
    if(JB_found==0){
        cout<<"JB_found=0\n";
        cout<<"stop"<<stop<<"\n";
        cout<<"cump"<<cump<<"\n";
        cout<<"sum"<<"\n";
        cout<<"urB.size()"<<urB.size()<<"\n";
        //cout<<"i "<<i<<endl;
        //cout<<"tmpP"<<tmpP<<"\n";
        system("pause");
    }
}


void SelectJunct_reverse(size_t &JA,size_t &JB)
{

    /*
    //Make array which can hold the possible junctions which can be picked for revrse rxn
    vector<size_t> av_rev_juncA, av_rev_juncB; // avaialable junctions for reverse rxn
    for (size_t i=0;i<reactedA.size();++i){ // all junctions in reacted array are eligible
        av_rev_juncA.push_back(reactedA[i]);
    }
    for (size_t i=0;i<urA.size();++i){// only those jucntions in unreated array are eligible which have at least one neighbor(ie at least one reacted node)
        if(neighA[urA[i]].size()>0){
            av_rev_juncA.push_back(urA[i]);}
    }
    for (size_t i=0;i<reactedB.size();++i){ // all junctions in reacted array are eligible
        av_rev_juncB.push_back(reactedB[i]);
    }
    for (size_t i=0;i<urB.size();++i){// only those jucntions in unreated array are eligible which have at least one neighbor(ie at least one reacted node)
        if(neighB[urB[i]].size()>0){
            av_rev_juncB.push_back(urB[i]);}
    }
    */

    // by randomly selecting bonds
   
    double rand_num=rand32()/4294967296.0;
    int bond_num=rand_num*all_bonds.size();
    JA=all_bonds[bond_num].JA;
    JB=all_bonds[bond_num].JB;
    all_bonds.erase(all_bonds.begin()+bond_num); // bond between JA and JB is broken now, so I can remove that bond from the list
    
}

void UpdateSum_forward(const size_t JA,const size_t JB,const size_t idxA,const size_t idxB)
{
    double probAB=0;
    for(size_t i=0;i<urA.size();++i) {
        tmpJA = urA[i];
        probAB = p[dist(JunctionDist[tmpJA][JB])] * (double)(fA - neighA[tmpJA].size());
        sum -= probAB;
        sumA[tmpJA] -= probAB;
    }
    for(size_t i=0;i<urB.size();++i) {
        tmpJB = urB[i];
        probAB = p[dist(JunctionDist[JA][tmpJB])] * (double)(fB - neighB[tmpJB].size());
        sum -= probAB;
        sumA[JA] -= probAB;
    }
    probAB = p[dist(JunctionDist[JA][JB])];
    sum += probAB;
    sumA[JA] += probAB;// avoid double counting
    // erase element for urA/urB
    // if size == fA-1 or fB-1, then all nodes have been reacted for the junction
    // neighbor list not updated here because it affects the calculation of original distance between junctions
    if(neighA[JA].size() == fA-1) {
        reactedA.push_back(JA);
        //cout<<"JA "<<JA<<" pushed to reactedA"<<endl;
        urA.erase(urA.begin()+idxA);
        //cout<<"JA "<<urA[idxA]<<" erased from unreactedA"<<endl;
        //if(JA==185) {cout<<"forward"<<endl;system("pause");}
    }
    if(neighB[JB].size() == fB-1) {
        reactedB.push_back(JB);
        urB.erase(urB.begin()+idxB);
    }
    
}



void UpdateSum_reverse(const size_t JA,const size_t JB)//,const size_t idxA,const size_t idxB)
{
    double sum_test=0;
    vector<double> sumA_test; 
    sumA_test = vector<double>(NRA,0);

    for(size_t i=0;i<urA.size();++i){
        size_t JA_loop=urA[i];
        double probAB;
        for(size_t j=0;j<urB.size();++j){
            size_t JB_loop=urB[j];
            probAB=p[dist(JunctionDist[JA_loop][JB_loop])]*(fA-(double)neighA[JA_loop].size())*(fB-(double)neighB[JB_loop].size());
            sum_test+=probAB;
            
            sumA_test[JA_loop]+=probAB;
            if(probAB<0){
                cout<<"probAB<0";
                system("pause");
            }
        }
    }
    sum=sum_test;
    sumA=sumA_test; // contains sum of probability of all unreacted junctions
    //sum_test=0; // reassigning value, but not redefining variable - because already defined
    //sumA_test = vector<double>(NRA,0);
}

void UpdateLoop_forward(const size_t JA,const size_t JB)
{
    for(size_t i=0;i<loop.size();++i) {
        if(dist(JunctionDist[JA][JB]) == i) {
        // loop[0] is for new branch, not loop  --> that is why ++i is done
        // for system with fA or fB =2, loop[n] is for n/2 order loop
        // for other system, loop[n] is for nth order loop
        // loop[n] should be zero for odd n
            loop[i]++;
            LoopSize[JA][JB]=i;
            break;
        }
    }
}


void UpdateLoop_reverse(const size_t JA,const size_t JB)
{
    //for(size_t i=0;i<loop.size();++i) {
        //if(mol[JA]==mol[JB+NRA]){
    //if(dist(JunctionDist[JA][JB]) == i) {
    // loop[0] is for new branch, not loop  --> that is why ++i is done
    // for system with fA or fB =2, loop[n] is for n/2 order loop
    // for other system, loop[n] is for nth order loop
    // loop[n] should be zero for odd n
    if(LoopSize[JA][JB]==1){
        loop[LoopSize[JA][JB]]--;// LoopSize[JA][JB]- is the loopsize in which JA and JB were involved before breaking
        //cout<<"Updated loop fraction of size"<<LoopSize[JA][JB]<<"\n";
    }
    LoopSize[JA][JB]=0;// because after bond breaking, there is no more loop between JA and JB- hence LoopSize=0

}

void CollectConnected(const size_t JA,const size_t JB,vector<size_t> &AcA,vector<size_t> &AcB,
vector<size_t> &BcA,vector<size_t> &BcB,vector<size_t> &dAcA,vector<size_t> &dBcB)
{   
    //cout<<"in collectconnected_forward"<<endl;
   // system("pause");
    int tmpDist;
    // collecting A junctions connected to JA/JB
    //cout<<"JunctionDist[278][188]"<<JunctionDist[278][188]<<endl;
    for(size_t i=0;i<NRA;++i) {
        tmpJA = i;
        // collecting A junction connected to JB
        if(JunctionDist[tmpJA][JB])
            AcB.push_back(tmpJA);
        // collecting A junction connected to JA
        tmpDist = getJunctionDistAA(tmpJA,JA);
        if(tmpDist >= 0) {
            AcA.push_back(tmpJA);
            dAcA.push_back((size_t)tmpDist);
        }
    }
    // collecting B junctions connected to JA/JB
    for(size_t i=0;i<NRB;++i) {
        tmpJB = i;
        // collecting B junction connected to JA
        if(JunctionDist[JA][tmpJB])
        BcA.push_back(tmpJB);
        // collecting B junction connected to JB
        tmpDist = getJunctionDistBB(tmpJB,JB);
        if(tmpDist >= 0) {
            BcB.push_back(tmpJB);
            dBcB.push_back((size_t)tmpDist);
        }
    }
    //AcB.push_back(JA); // not actually needed 
   // BcA.push_back(JB);
    /*cout<<"BcB[j]:";
    for(size_t j=0;j<BcB.size();++j) cout<<BcB[j]<<" ";
        cout<<"BcA[j]:";
    for(size_t j=0;j<BcA.size();++j) cout<<BcA[j]<<" ";
        cout<<endl;
        cout<<"AcA[j]:";
    for(size_t j=0;j<AcA.size();++j) cout<<AcA[j]<<" ";
        cout<<endl;
        cout<<"AcB[j]:";
    for(size_t j=0;j<AcB.size();++j) cout<<AcB[j]<<" ";
        cout<<endl;
    cout<<"JunctionDist[JA][JB] before update"<<JunctionDist[JA][JB]<<endl;*/
    
    /*if((AcB.size()==0 && BcA.size()!=0)|| (BcA.size()==0 && AcB.size()!=0)){
        cout<<"Problem in collect_connected-forward- AcB-BcA"<<endl;
        
        system("pause");
    //system("pause");
    }
    if((BcB.size()>0 && AcB.size()==0)){
        cout<<"Problem in collect_connected-forward- BcB-AcB"<<endl;
        
        system("pause");
    //system("pause");
    }*/
    
    //system("pause");
}




void CollectConnected_reverse(const size_t JA,const size_t JB,vector<size_t> &AcA,vector<size_t> &AcB,
vector<size_t> &BcA,vector<size_t> &BcB,vector<size_t> &dAcA,vector<size_t> &dBcB)
{
    /*if(JA==328 && JB==104){
        cout<<"BcB[j]:";
        for(size_t j=0;j<BcB.size();++j) cout<<BcB[j]<<" ";
        cout<<"BcA[j]:";
        for(size_t j=0;j<BcA.size();++j) cout<<BcA[j]<<" ";
        cout<<endl;
        cout<<"AcA[j]:";
        for(size_t j=0;j<AcA.size();++j) cout<<AcA[j]<<" ";
        cout<<endl;
        cout<<"AcB[j]:";
        for(size_t j=0;j<AcB.size();++j) cout<<AcB[j]<<" ";
        cout<<endl;
        system("pause");
    }*/
    int tmpDist;
    // collecting A junctions connected to JA/JB
    for(size_t i=0;i<NRA;++i) {
        tmpJA = i;
        // collecting A junction connected to JB
        if(JunctionDist[tmpJA][JB]>0 && JunctionDist[tmpJA][JB]!=MAX)// && tmpJA!=JA)
            AcB.push_back(tmpJA);
            // collecting A junction connected to JA
        tmpDist = getJunctionDistAA(tmpJA,JA);
        if(tmpDist >= 0) {// will never be equal to zero except for self
            AcA.push_back(tmpJA);
            dAcA.push_back((size_t)tmpDist);
        }
    }
    // collecting B junctions connected to JA/JB
    for(size_t i=0;i<NRB;++i) {
        tmpJB = i;
        // collecting B junction connected to JA
        if(JunctionDist[JA][tmpJB]>0 && JunctionDist[tmpJA][JB]!=MAX)// && tmpJB!=JB)
            BcA.push_back(tmpJB);
            // collecting B junction connected to JB
        tmpDist = getJunctionDistBB(tmpJB,JB);
        if(tmpDist >= 0) {
            BcB.push_back(tmpJB);
            dBcB.push_back((size_t)tmpDist);
        }
    }
    //cout<<"In CollectConnected_reverse\n";
    /*for (size_t k=0; k< neighB[JB].size(); ++k){
    if (std::count(AcB.begin(), AcB.end(),neighB[JB][k])) {
            //std::cout << "neighB[JB][k]"<<neighB[JB][k]<<" found in AcB\n";
            
        }
    else{
        std::cout << "neighB[JB][k] NOT found in AcB\n";
        cout<<"JunctionDist[neighB[JB][k]][JB]"<<JunctionDist[neighB[JB][k]][JB]<<"\n";
        cout<<"k"<<k<<"\n";
        system("pause");
    }
    }*/


    /*for (size_t k=0; k< neighA[JA].size(); ++k){
    if (std::count(BcA.begin(), BcA.end(),neighA[JA][k])) {
            //std::cout << "neighA[JA][k] found in BcA\n";
            
        }
    else{
        std::cout << "neighA[JA][k] NOT found in BcA\n";
        cout<<"JunctionDist[JA][neighA[JA][k]]"<<JunctionDist[JA][neighA[JA][k]]<<"\n";
        cout<<"k"<<k<<"\n";
        system("pause");
    }
    }*/
     /*cout<<"BcB[j]:";
    for(size_t j=0;j<BcB.size();++j) cout<<BcB[j]<<" ";
        cout<<"BcA[j]:";
    for(size_t j=0;j<BcA.size();++j) cout<<BcA[j]<<" ";
        cout<<endl;
        cout<<"AcA[j]:";
    for(size_t j=0;j<AcA.size();++j) cout<<AcA[j]<<" ";
        cout<<endl;
        cout<<"AcB[j]:";
    for(size_t j=0;j<AcB.size();++j) cout<<AcB[j]<<" ";
        cout<<endl;
    cout<<"JunctionDist[JA][JB] in collectconnected_reverse: "<<JunctionDist[JA][JB]<<endl;*/

   // if((AcB.size()==0 && BcA.size()!=0)|| (BcA.size()==0 && AcB.size()!=0)){
   //    cout<<"Problem in collect_connected-reverse: AcB-BcA"<<endl;
    //system("pause");}
    /*if((BcB.size()>1 && AcB.size()==0)){
        cout<<"Problem in collect_connected-reverse- BcB-AcB"<<endl;
        
        system("pause");
    //system("pause");
    }
    if((AcA.size()>1 && BcA.size()==0)){ // AcA size will always be at leats 1 because it includes JA itself
        cout<<"Problem in collect_connected-reverse- AcA-BcA"<<endl;
        
        system("pause");
    //system("pause");
    }
    */

}

void UpdateJuncDist_forward(const size_t JA,const size_t JB,const vector<size_t> &BcA,const vector<size_t> &AcB)
{
   // cout<<"JA_for: "<<JA<<", "<<"JB_for: "<<JB<<endl;

    size_t olddist,newdist;
    double probAB;
    for(size_t i=0;i<BcA.size();++i) {
        tmpJB = BcA[i];
        for(size_t j=0;j<AcB.size();++j) {
            tmpJA = AcB[j];
            olddist = JunctionDist[tmpJA][tmpJB];
            newdist = JunctionDist[tmpJA][JB] + JunctionDist[JA][tmpJB] + 1;
            if(olddist == 0) { // i.e. previously unconnected
            if(newdist > MAX) newdist = MAX;
            JunctionDist[tmpJA][tmpJB] = newdist;
            //if distance is updated, sum of relative probabilities must also be updated
            probAB = (p[dist(newdist)] - p[dist(olddist)]) * (double)(fA-neighA[tmpJA].size())*(double)(fB-neighB[tmpJB].size());
            if( (neighA[tmpJA].size()<fA) && (neighB[tmpJB].size()<fB) ) {
                sum += probAB;
                sumA[tmpJA] += probAB;
            }
            continue;
            }
            if(newdist < olddist) {
                JunctionDist[tmpJA][tmpJB] = newdist;
                //if distance is updated, sum of relative probabilities must also be updated
                probAB = (p[dist(newdist)] - p[dist(olddist)]) * (double)(fA-neighA[tmpJA].size())*(double)(fB-neighB[tmpJB].size());
                if( (neighA[tmpJA].size()<fA) && (neighB[tmpJB].size()<fB) ) {
                    sum += probAB;
                    sumA[tmpJA] += probAB;
                }
                continue;
            }
        }
    }
}

void UpdateJuncDist_reverse(const size_t JA,const size_t JB)//,const vector<size_t> &BcA,const vector<size_t> &AcB, const vector<vector<size_t> >&neighA,const vector<vector<size_t> >&neighB)
{

    // IMPORTANT comment for debugging
    /*for(size_t i=0;i<NRA;++i){
        cout<<"i"<<i<<endl;
        for (size_t j=0;j<neighA[i].size();++j){
            cout<<neighA[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<"JA_rev: "<<JA<<", "<<"JB_rev: "<<JB<<endl;
*/


   /* cout<<"JunctionDist before reverse update:"<<endl;
    for (int i=0;i<NRA;++i){
        for(int j=0;j<NRB;++j){
            cout<<JunctionDist[i][j]<<" ";
        }
        cout<<endl;
    }
    system("pause");*/
    
    //cout<<"inside update junction distance_reverse\n";
    //system("pause");
    size_t olddist,newdist;
    //olddist=NRA+NRB; // just for initializing
    double probAB;
    // since this is for reverse reaction: we have to update the JunctionDist matrix using floyd warshall algo
    // since this also involves updating the sum and probAB- I might just implement floyd-warshall in this step directly instead of making a whole new function
    // so- the main difference between reverse sum update(floyd warshall) and forward sum update is that 
    // for forward- I can just check the distance JunctionDist[tmpJA][JB] + JunctionDist[JA][tmpJB] + 1 for the junctions JA anf JB only - which are reacting in that step
    // but for reverse- I will have to check this distance for all junctions and not only JA and JB- because of the inherent loss in shortest distance info when a bond is broken
    // so this adds another loop where JA and JB are being varied
    // this is basically floyd -warshall!!
    vector <size_t> junctions_clusterA; // lists all the jucntions present in the cluster/molecule to which A and dB belonged to before unreacting 
    vector <size_t> junctions_clusterB;
    size_t molecule=mol[JA];// this molecule will be the same for A and B
    //cout<<"moleculeA"<<mol[JA]<<"\n";
    //cout<<"moleculeB"<<mol[JB+NRA]<<"\n";
    //if(mol[JA]!=mol[NRA+JB]){cout<<"Problem in update junction reverse\n"; system("pause");}
    int count=0; // counts the number of tiems I encounter JA_loop=Ja and JB_loop=JB
    //cout<<"molecule JA"<<mol[JA]<<"\n";
    //cout<<"molecule JB"<<mol[JB+NRA]<<"\n";
    //cout<<"junctions_cluster\n";
    //cout<<"JA inside updatejunction: "<<JA<<"\n";
    //cout<<"JB inside updatejunction: "<<JB<<"\n";

    for (size_t i=0;i<NRA; ++i){
        if(mol[i]==molecule) {junctions_clusterA.push_back(i);}// cout<<i<<" ";} // find all A junctions contained in the cluster to which JA belongs
    } // equivalent to AcA- before reverse update
    for (size_t i=0;i<NRB; ++i){
        if(mol[i+NRA]==molecule) {junctions_clusterB.push_back(i); }//cout<<i<<" ";} // find all A junctions contained in the cluster to which JA belongs
    }
    size_t loop_exist=0;
    for (size_t i=0;i<neighA[JA].size();++i){
        if(neighA[JA][i]==JB) {//JunctionDist[JA][JB]=1;// since they are neighbours, junctions are NOT being temporariliy disconnected!!
        //cout<<"JA and JB are neighbours, hence not temporarily disconnected!!\n";
        //cout<<"LOOP- neighbors even after disconnected\n";
        loop_exist=1;
        break;
        }
    }
    if(loop_exist==0){ // no loop breaking- hence will have to update junction distances
        
        vector<vector<unsigned short> > JunctionDist_copy;
        JunctionDist_copy=JunctionDist;
        JunctionDist[JA][JB]=0;// temporarily disconnected jucntions JA and JB
        //cout<<"H1\n";
        for(size_t i=0;i<junctions_clusterA.size();++i){
            size_t JA_loop=junctions_clusterA[i];
            for(size_t j=0;j<junctions_clusterB.size();++j){
                size_t JB_loop=junctions_clusterB[j];
                //cout<<"JB_loop in cluster:"<<JB_loop<<"\n";
            //for(size_t j=0;j<junctions_clusterB.size();++j){
            //  size_t JB_loop=junctions_clusterB[j];
                if((JunctionDist[JA_loop][JB_loop]>1)){//&& (JunctionDist[JA_loop][JB]+JunctionDist[JA][JB_loop]+1!=1)){//there exists an indirect connection through [JA][JB] not direct connection!!- similar to rui's code
                   // cout<<"i:"<<i<<endl;
                    //cout<<"JA_loop: "<<JA_loop<<" "<<"JB_loop: "<<JB_loop<<endl;
                    //cout<<"neighA[JA_loop].size()"<<neighA[JA_loop].size()<<endl;
                   // cout<<"JunctionDist[JA_loop][JB_loop]"<<JunctionDist[JA_loop][JB_loop]<<endl;
                    //system("pause");
                    JunctionDist[JA_loop][JB_loop]=0;
                    
                }
            }
        }
    vector<vector<unsigned short> > JunctionDist_FW=JunctionDist;
    //cout<<"junctions_clusterA.size()"<<junctions_clusterA.size()<<endl;
    //cout<<"junctions_clusterB.size()"<<junctions_clusterB.size()<<endl;

    /*cout<<"JunctionDist before reverse update but after modification:"<<endl;
    for (int i=0;i<junctions_clusterA.size();++i){
        for(int j=0;j<junctions_clusterB.size();++j){
            cout<<JunctionDist[junctions_clusterA[i]][junctions_clusterB[j]]<<" ";
        }
        cout<<endl;
    }
    system("pause");*/

    /*
    cout<<"junctions_clusterB:"<<endl;
    for (int i=0;i<junctions_clusterB.size();++i){
        cout<<junctions_clusterB[i]<<" ";
    }
    cout<<endl;
    cout<<"junctions_clusterA:"<<endl;
    for (int i=0;i<junctions_clusterA.size();++i){
        cout<<junctions_clusterA[i]<<" ";
    }
    cout<<endl;
    */
    //system("pause");

// n-BFS starts here!!
vector<vector< int> > case_number_BFS;
vector<vector< int> > case_number_FW;
for(size_t i=0;i<NRA;++i){
        case_number_BFS.push_back(vector<int>(NRB,-1));
        //JunctionDist_FW.push_back(vector<unsigned short>(NRB,0));

        //vector<vector<unsigned short> > JunctionDist_FW;

        case_number_FW.push_back(vector<int>(NRB,-1));
        
    }



// for BFS- I only need the connectivity info about crosslinks, and not the chains, because the chains act as the edges of the graph, which is already present in connectivity info of crosslinks
    if(junctions_clusterA.size()>1){    
        for(size_t A_idx=0;A_idx<junctions_clusterA.size();++A_idx){ // s- source node # do BFS for each source node s
            //cout<<"NEW A_idx:"<<A_idx<<endl;
           // cout<<endl;
            //cout<<endl;
           // system("pause");
            size_t Anode=junctions_clusterA[A_idx];
            
        // single BFS run starts here
        // this will give me the shortest distance between two crosslinkers(B nodes)
        // and then I will have to convert that into the shortest distance between 
        //cout<<"neighA[Anode].size()"<<neighA[Anode].size()<<endl;
        if(neighA[Anode].size()>0){
        // the BFS algo assumes that the graph is connected
        // so, starting from source node s- the algorithm will only update the shortest distances within the cluster of connected nodes
            size_t s=neighA[Anode][0]; // source B node
            
            //cout<<"s"<<s<<endl;
            //cout<<"Anode"<<Anode<<endl;
            //system("pause");
            size_t visited1[NRB]={0};
            int d1[NRB]; // distance - this is slightly higher on the memory- because not all the crosslinks are present in this cluster, but I am using the full array just to keep the indices consistent
            int d[NRB]; 
            for (size_t i=0;i<NRB;++i){
                d1[i]=-1;
                d[i]=-1;
            }
            visited1[s]=1;
            d1[s]=0;
            queue<size_t> q1;
            q1.push(s);
            while(!q1.empty()){
                int v1=q1.front();
                /*cout<<"v1:"<<v1<<endl;
                cout<<"G[v1].size():"<<G[v1].size()<<endl;
                //cout<<"G[1].size():"<<G[1].size()<<endl;
                cout<<"G.size():"<<G.size()<<endl;
                cout<<"junctions_clusterB.size():"<<junctions_clusterB.size()<<endl;
                cout<<"G[v1]:"<<endl;

                for (int i=0;i<G[v1].size();++i){
                    cout<<G[v1][i]<<" ";
                }
                cout<<endl;*/
                
                //system("pause");
                //cout<<endl;
                q1.pop();
                for(size_t u1: G[v1]){
                    if(!visited1[u1]){
                        visited1[u1]=1;
                        q1.push(u1);
                        d1[u1]=d1[v1]+1;
                        //cout<<"u1="<<u1<<endl;
                        //cout<<"d1[u1]="<<d1[u1]<<endl;
                    }
                }
            }
int d2[NRB]; // distance - this is slightly higher on the memory- because not all the crosslinks are present in this cluster, but I am using the full array just to keep the indices consistent
            if(neighA[Anode].size()>1){
                
                size_t s2=neighA[Anode][1]; // hardcoded for bifunctional polymer  
                size_t visited2[NRB]={0};
                
                for (size_t i=0;i<NRB;++i){
                    d2[i]=-1;
                }
                visited2[s2]=1;
                d2[s2]=0;
                queue<size_t> q2;
                q2.push(s2);
                while(!q2.empty()){
                    int v2=q2.front();
                    q2.pop();
                    for(size_t u2: G[v2]){
                        if(u2==neighA[Anode][0]){ // optimizing such that the second one doesn't go over the same path again, but with larger distance because of overlap with chain
                            continue; // continue here because we need to look at other paths!!!!!

                        }
                        else if(!visited2[u2]){
                            visited2[u2]=1;
                            q2.push(u2);
                            d2[u2]=d2[v2]+1;
                            //cout<<"d2[u2]"<<d2[u2]<<endl;
                        }
                    }
                }

                for(size_t i=0;i<NRB;++i){
                    if(d1[i]==-1 && d2[i]!=-1){
                         d[i]=d2[i];
                          

                    }
                    else if (d2[i]==-1 && d1[i]!=-1) d[i]=d1[i];
                    else if(d1[i]!=0 && d2[i]!=0 && d1[i]!=-1 && d2[i]!=-1){
                        d[i]=min(d1[i],d2[i]);
                        /*if(d[i]>1000){
                            cout<<"big value of d assigned:"<<d[i]<<endl;
                            system("pause");
                    }*/
                        
                    }
                    else{
                       // cout<<"i"<<i<<endl;
                       // cout<<"d1[i]"<<d1[i]<<endl;
                        //cout<<"d2[i]"<<d2[i]<<endl;
                        //cout<<"d1 or d2 is zero or -1"<<endl;
                        if(d1[i]==0 && d2[i]==0){
                            d[i]=d1[i];
                           
                        }
                        else{
                            //cout<<"d value not assigned!-PROBLEM MAYBE"<<endl;
                            //system("pause");
                        }
                       
                    }
                    
                }
            }
            else{
                for(size_t i=0;i<NRB;++i){
                    d[i]=d1[i];
                }
            }
            
            //else{
                for (int i=0;i<junctions_clusterB.size();++i){
                    size_t temp=junctions_clusterB[i];
                    //if(){}
                   /* if(Anode==2 && temp==0){
                            cout<<"Anode==2 && temp==0"<<endl;
                            cout<<"JunctionDist[Anode][temp]"<<JunctionDist[Anode][temp]<<endl;
                            cout<<"d1[temp]"<<d1[temp]<<endl;
                            cout<<"d2[temp]"<<d2[temp]<<endl;
                            cout<<"d[temp]"<<d[temp]<<endl;
                            system("pause");
                        }*/
                    //cout<<"d[temp]: "<<d[temp]<<endl;
                    if(JunctionDist[Anode][temp]!=1){ // whatever was directly connected stays directly connected- no need to change
                        /*if(d[temp]>1000){
                                    cout<<"Anode"<<Anode<<endl;
                                    cout<<"temp"<<temp<<endl;
                                    cout<<"big value of d"<<endl;
                                    cout<<"d1[temp]"<<d1[temp]<<endl;
                                    cout<<"d2[temp]"<<d2[temp]<<endl;
                                    cout<<"d[temp]"<<d[temp]<<endl;
                                    system("pause");
                                }*/
                    if(d[temp]==0){
                        case_number_BFS[Anode][temp]=1;
                        //cout<<"JunctionDist[Anode][temp]"<<JunctionDist[Anode][temp]<<endl;

                        JunctionDist[Anode][temp]=d[temp]; // completely disconnected
//  cout<<"case number=1"<<endl;
// cout<<"Anode"<<Anode<<endl;
                        //cout<<"temp"<<temp<<endl;
                            // cout<<"d[temp]: "<<d[temp]<<endl;
                       // temp_JuncDist_array[tmpcnt][0]=temp; // the corresponding index
                        //temp_JuncDist_array[tmpcnt][1]=d[temp]; // the corresponsing value
                        //tmpcnt++;
                        //cout<<"Completely disconnected now"<<endl;
                        //cout<<"Anode"<<Anode<<" temp"<<temp<<endl;
                        //system("pause");    
                    }   
                        //cout<<"d=0- NOT ALLOWED in BFS"<<endl;
                        //system("pause");
                        /*else{
                            // happens for primary loops!!
                            cout<<"PROBLEM!! trying to update distance of directly connected nodes"<<endl;
                            cout<<"Anode"<<Anode<<" temp"<<temp<<endl;
                            system("pause");
                        }*/
                        //JunctionDist[Anode][temp]=d1[temp];
                    
                
                    else if (d[temp]!=-1 ){
                        case_number_BFS[Anode][temp]=2;
                       // if(JunctionDist[Anode][temp]!=1){// or JunctionDist[Anode][temp]!=0)

                            JunctionDist[Anode][temp]=2*d[temp]+1; // 1 is  added because have to consider the extra chain
                            
                            //temp_JuncDist_array[tmpcnt][0]=temp; // the corresponding index
                            //temp_JuncDist_array[tmpcnt][1]=2*d[temp]+1; // the corresponsing value
                            //tmpcnt++;
                            // cout<<"case number=2"<<endl;
// cout<<"Anode"<<Anode<<endl;
                        //cout<<"temp"<<temp<<endl;
                            //   cout<<"d[temp]: "<<d[temp]<<endl;
                            //cout<<"Anode:"<<Anode<<endl;
                            //cout<<"temp:"<<temp<<endl;
                       // }
                        /*else{
                            cout<<"PROBLEM!! trying to update distance of directly connected nodes- maybe loop"<<endl;
                            cout<<"Anode"<<Anode<<" temp"<<temp<<endl;
                            system("pause");
                        }*/
                    }
                    else{
                        case_number_BFS[Anode][temp]=3;
                        //cout<<"distances not changed!! why??"<<endl;
                        //system("pause");
                    }
                }
            }
            //}
        }

        }
    }

/*
// Floyd-Warshall starts here!
//size_t case_number_FW=-1;

for(size_t ite=0;ite<all_bonds.size();++ite){ // dont put JA or JB in loop- because they are global variables and value will get changed without you even knowing what is happening!!!
            size_t JA_k=all_bonds[ite].JA; // JA_loop (s) are the A jucntions which have to be looped over in order to update the shortest distance
            size_t JB_k=all_bonds[ite].JB;

            for(size_t i=0;i<junctions_clusterA.size();++i) {
                //tmpJB = BcA[i];
                tmpJA=junctions_clusterA[i];
                for(size_t j=0;j<junctions_clusterB.size();++j) {
                    tmpJB = junctions_clusterB[j];
                    olddist = JunctionDist_FW[tmpJA][tmpJB]; 
                    if(JunctionDist_copy[tmpJA][tmpJB]!=0 && (tmpJA!=JA_k) && (tmpJB!=JB_k)){//if something was disconnected before rev rxn- it stays disconnected
                        if(JunctionDist_FW[tmpJA][JB_k]!=0 && JunctionDist_FW[JA_k][tmpJB]!=0 && JunctionDist_FW[JA_k][JB_k]!=0) // checking for infinity cases
                            newdist = max(JunctionDist_FW[tmpJA][JB_k] + JunctionDist_FW[JA_k][tmpJB]+JunctionDist_FW[JA_k][JB_k],0);// JunctionDist[JA_loop][JB_loop]: this additional term was =1 in case of forward- here i dont know what this is!!
                        else newdist=0;
                        // IMP: since this is breakage- olddist=0 means temporarily disconnected onyl
                        // there can be NO other scanario where olddist=0
                        if(JunctionDist_FW[tmpJA][tmpJB] == 0 && newdist!=0) { // this has to be JunctionDist[tmpJA][tmpJB] and not JunctionDist_copy[tmpJA][tmpJB]
                            // this means- junctions were temporarily unconnected due to breakage
                        // but in this case- sum still needs to be updated because physically- newdist is less than olddist (which was essentially infinity in physical sense, but represented as 0- so the math doesn't really work so well there!!!)
                            if(newdist==1) {newdist=0;
                            cout<<"in newdist=1 case\n";
                            system("pause");
                            } //olddist=0 means temp disconnected- and newdist=1 on top of that means that the actual connection was only a direct one, and no other loop type connection exists-
                        // in this case- they will be actually  disconnecetd after breakage
                            if(newdist > MAX) newdist = 0;
                            JunctionDist_FW[tmpJA][tmpJB] = newdist;
                            case_number_FW[tmpJA][tmpJB]=1;
                        }
                        else if(newdist==0 && olddist>0){ // this case would mean that there is some path by which the two junctions become unconnected, even thoguh they were connected by some other path
                        // in this case- this new dist is NOT VALID!!-YESS EXACTLY- thats why the the sum update step also needs to be modified!
                        //basically- in this case- we dont need the sum update step because distance was never updated in the first place!
                        // this extra loop is required because I am denoting unconnectedness by 0 and not by infinity (practically some max number!!) 
                        
                        // this means that the dist is actually oldist only, just that for this JA_k and JB_k- the dist turns out to be 0, but that should not be trusted!
                          //cout<<"in newdist=0 case\n";
                           // system("pause");
                        // JunctionDist[tmpJA][tmpJB] = olddist;
                            //cout<<"into this extra if condition\n";
                            //cout<<"JunctionDist["<<tmpJA<<"]["<<tmpJB<<"]="<<JunctionDist[tmpJA][tmpJB]<<"\n";

                        //cout<<"JunctionDist["<<tmpJA<<"]["<<tmpJB<<"="<<newdist<<"\n";
                        //if distance is updated, sum of relative probabilities must also be updated
                        //probAB = (p[dist(newdist)] - p[dist(olddist)]) * (double)(fA-neighA[tmpJA].size())*(double)(fB-neighB[tmpJB].size());

                        }
                        else if(newdist < olddist) {// both non-zero
                            JunctionDist_FW[tmpJA][tmpJB] = newdist;
                            case_number_FW[tmpJA][tmpJB]=2;
                        }
                    }
                    //} if condition bracket commented
                }
                }

        }

for (int i=0;i<NRA;++i){
        for(int j=0;j<NRB;++j){
            if(JunctionDist[i][j]!=JunctionDist_FW[i][j]){
                cout<<"ERROR:: JunctionDist[i][j]!=JunctionDist_FW[i][j] !!"<<endl;
                cout<<"i: "<<i<<endl;
                cout<<"j: "<<j<<endl;
                cout<<"JunctionDist_FW[i][j]: "<<JunctionDist_FW[i][j]<<endl;
                cout<<"JunctionDist_BFS[i][j]: "<<JunctionDist[i][j]<<endl;
                cout<<"JunctionDist_copy[i][j]: "<<JunctionDist_copy[i][j]<<endl;
                cout<<"JA: "<<JA<<endl;
                cout<<"JB: "<<JB<<endl;
                cout<<"case_number_FW["<<i<<"]["<<j<<"]: "<<case_number_FW[i][j]<<endl;
                cout<<"case_number_BFS["<<i<<"]["<<j<<"]: "<<case_number_BFS[i][j]<<endl;

                cout<<"neighA[i]:"<<endl;
                for (int k=0;k<neighA[i].size();++k){
                    cout<<neighA[i][k]<<" ";
                }
                cout<<endl;

                cout<<"neighB[j]:"<<endl;
                for (int k=0;k<neighB[j].size();++k){
                    cout<<neighB[j][k]<<" ";
                }
                cout<<endl;
                cout<<endl;
                cout<<"JA: "<<JA<<endl;
                cout<<"JB: "<<JB<<endl;
                cout<<"neighA[JA]:"<<endl;
                for (int k=0;k<neighA[JA].size();++k){
                    cout<<neighA[JA][k]<<" ";
                }
                cout<<endl;

                cout<<"neighB[JB]:"<<endl;
                for (int k=0;k<neighB[JB].size();++k){
                    cout<<neighB[JB][k]<<" ";
                }
                cout<<endl;

                cout<<"G[JB]:"<<endl;
                for (int k=0;k<G[JB].size();++k){
                    cout<<G[JB][k]<<" ";
                }

                cout<<endl;
                cout<<"G[j]:"<<endl;
                for (int k=0;k<G[j].size();++k){
                    cout<<G[j][k]<<" ";
                }

                //cout<<"neighA[j].size()= ";
                //cout<<neighA[i].size()<<endl;
                //cout<<"neighB[i].size()= ";
                //cout<<neighB[j].size()<<endl;
                //cout<<"neighA.size()= "<<neighA.size()<<endl;
                //cout<<"neighB.size()= "<<neighB.size()<<endl;

                for(int k=0;k<neighA.size();++k){
                    cout<<"k="<<k<<" neighA[k].size= "<<neighA[k].size()<<endl;
                }
                for(int k=0;k<neighB.size();++k){
                    cout<<"k="<<k<<" neighB[k].size= "<<neighB[k].size()<<endl;
                }

           
                
                cout<<endl;

                

                system("pause");
            }
        }
        //cout<<endl;
    }

*/

    }
    /*cout<<"JunctionDist_reverse:"<<endl;
    for (int i=0;i<NRA;++i){
        for(int j=0;j<NRB;++j){
            cout<<JunctionDist[i][j]<<" ";
        }
        cout<<endl;
    }*/
    //cout<< (double)(clock()-c0)/CLOCKS_PER_SEC<<endl;
    //system("pause");
}



void UpdateJuncDist_common(const size_t JA,const size_t JB,const vector<size_t> &AcA,const vector<size_t> &BcB,
const vector<size_t> &dAcA,const vector<size_t> &dBcB)
{
    size_t olddist,newdist;
    double probAB;
    for(size_t i=0;i<AcA.size();++i) {
        tmpJA = AcA[i];
        for(size_t j=0;j<BcB.size();++j) {
            tmpJB = BcB[j];
            olddist = JunctionDist[tmpJA][tmpJB];
            newdist = dAcA[i] + dBcB[j] + 1;
            if(olddist == 0) { // i.e. previously unconnected
                if(newdist > MAX) newdist = MAX;
                JunctionDist[tmpJA][tmpJB] = newdist;
                //if distance is updated, sum of relative probabilities must also be updated
                probAB = (p[dist(newdist)] - p[dist(olddist)]) * (double)(fA-neighA[tmpJA].size())*(double)(fB-neighB[tmpJB].size());
                if( (neighA[tmpJA].size()<fA) && (neighB[tmpJB].size()<fB) ) {
                    sum += probAB;
                    sumA[tmpJA] += probAB;
                }
                continue;
            }
            if(newdist < olddist) {
                JunctionDist[tmpJA][tmpJB] = newdist;
                //if distance is updated, sum of relative probabilities must also be updated
                probAB = (p[dist(newdist)] - p[dist(olddist)]) * (double)(fA-neighA[tmpJA].size())*(double)(fB-neighB[tmpJB].size());
                if( (neighA[tmpJA].size()<fA) && (neighB[tmpJB].size()<fB) ) {
                    sum += probAB;
                    sumA[tmpJA] += probAB;
                }
                continue;
            }
        }
    }
}

void UpdateJuncDist_common_reverse(const size_t JA,const size_t JB,const vector<size_t> &AcA,const vector<size_t> &BcB,
const vector<size_t> &dAcA,const vector<size_t> &dBcB)
{
    size_t olddist,newdist;
    double probAB;
    for(size_t i=0;i<AcA.size();++i) {
        tmpJA = AcA[i];
        for(size_t j=0;j<BcB.size();++j) {
            tmpJB = BcB[j];
            olddist = JunctionDist[tmpJA][tmpJB];
            newdist = dAcA[i] + dBcB[j] + 1;
            if(olddist == 0) { // i.e. previously unconnected
                if(newdist > MAX) newdist = MAX;
                //////JunctionDist[tmpJA][tmpJB] = newdist;
                //if distance is updated, sum of relative probabilities must also be updated
                probAB = (p[dist(newdist)] - p[dist(olddist)]) * (double)(fA-neighA[tmpJA].size())*(double)(fB-neighB[tmpJB].size());
                if( (neighA[tmpJA].size()<fA) && (neighB[tmpJB].size()<fB) ) {
                ////sum += probAB;
                ////sumA[tmpJA] += probAB;
                }
                continue;
            }
            if(newdist < olddist) {
            //////JunctionDist[tmpJA][tmpJB] = newdist;
            //if distance is updated, sum of relative probabilities must also be updated
                probAB = (p[dist(newdist)] - p[dist(olddist)]) * (double)(fA-neighA[tmpJA].size())*(double)(fB-neighB[tmpJB].size());
                if( (neighA[tmpJA].size()<fA) && (neighB[tmpJB].size()<fB) ) {
            ////sum += probAB;
            ////sumA[tmpJA] += probAB;
                }
            continue;
            }
        }
    }
}

void UpdateMol_forward(const size_t JA,const size_t JB,const vector<size_t> &AcA,const vector<size_t> &AcB,
const vector<size_t> &BcA,const vector<size_t> &BcB)
{
    size_t mA,mB,newMol,delMol;
    mA = mol[JA];// mol is molecule index for each junction(monomer) ), so mA is molecule index of JA, and mB is molecule index of JB (since this is B, NRA+JB)
    mB = mol[NRA+JB];
    //cout<<"mA "<<mA<<", mB "<<mB<<"\n";

    // new molecule number should be the smaller of the two molecule numbers
    if(mA != mB) { // if mA==mB then this is intramolecular rxn, no change is needed
    newMol = (mA<mB)?mA:mB;// min(mA,mB)
    delMol = (mA<mB)?mB:mA;// max(mA,mB)
    //cout<<"newMol "<<newMol<<"\n";
    //cout<<"delMol "<<delMol<<"\n";
    for(size_t j=0;j<AcA.size();++j) mol[AcA[j]] = newMol;
    for(size_t j=0;j<BcA.size();++j) mol[NRA+BcA[j]] = newMol;
    for(size_t j=0;j<AcB.size();++j) mol[AcB[j]] = newMol;
    for(size_t j=0;j<BcB.size();++j) mol[NRA+BcB[j]] = newMol;
    //molSize[newMol] = molSize[mA] + molSize[mB];
    //molSize[delMol] = 0;
    mol[JA]=newMol;
    mol[JB+NRA]=newMol;
    //if(molSize[newMol]> molSize[largestMol]) largestMol = newMol;

    }
    //if(JA==328 && mA==606) {cout<<"forward";system("pause");}

}


void UpdateMol_reverse(const size_t JA,const size_t JB,const vector<size_t> &AcA,const vector<size_t> &AcB,
const vector<size_t> &BcA,const vector<size_t> &BcB){
size_t mA,mB,newMol,oldMol;
mA = mol[JA];// mol is molecule index for each junction(monomer) ), so mA is molecule index of JA, and mB is molecule index of JB (since this is B, NRA+JB)
mB = mol[NRA+JB];
/*cout<<"JA"<<JA<<"JB"<<JB<<"\n";
cout<<"mA"<<mA<<"mB"<<mB<<"\n";
cout<<"JunctionDist[JA][JB]"<<JunctionDist[JA][JB]<<"\n";
cout<<"all_bonds.size()"<<all_bonds.size()<<"\n";
//system("pause");
cout<<"all_bonds:"<<"\n";

for (int i=0; i<all_bonds.size(); ++i){
    cout<<"JA"<<all_bonds[i].JA<<"JB"<<all_bonds[i].JB<<"\n";
}

*/
/*cout<<"mA"<<mA<<"mB"<<mB<<"\n";
cout<<"molSize before:\n";
for (size_t i=0;i<molSize.size(); i++){
    cout<<molSize[i]<<" ";
}
cout<<"\n";*/
// new molecule number should be the smaller of the two molecule numbers
//if(mA==0 && mB!=0){cout<<"mA=0 && mB!=0, so There might be problem right here!\n";system("pause");}
//cout<<"Inside UpdateMol_reverse\n";
//cout<<"JunctionDist[JA][JB]"<<JunctionDist[JA][JB]<<"\n";
size_t JA_JB_neigh=0;
for (int i=0;i<neighA[JA].size();++i){
    if(neighA[JA][i]==JB){
        JA_JB_neigh=1;
       // cout<<"JA and JB are neighbors\n";
        break;
    }
}
if(mA!=mB){cout<<"Error! reverse reactions species belong to different clusters\n"; 

cout<<"mA"<<mA<<"mB"<<mB<<"\n";
cout<<"JA"<<JA<<"JB"<<JB<<"\n";
cout<<"JunctionDist[JA][JB]"<<JunctionDist[JA][JB]<<"\n";
cout<<"\n";
system("pause");
}
    else if(JA_JB_neigh==0){ // if mA==mB then this is intramolecular rxn, no change is needed
        if(JunctionDist[JA][JB]==0){ // implies that JA and JB have been disconnected now, otherwise- no change because cluster is not changing
            oldMol=mol[JA];
            //cout<<"oldMol: "<<oldMol<<":\n";
            newMol=*max_element(mol.begin(), mol.end())+1;
            if(newMol>mol.size()){
                //cout<<"newMol>mol.size()\n";
                //system("pause");
            }
            //cout<<"newMol: "<<newMol<<":\n";

            mol[NRA+JB]=newMol;
           // cout<<"AcA[j]:\n";

            // THIS IS CORRECT!!!
            for(size_t j=0;j<AcA.size();++j) {mol[AcA[j]] = oldMol; //cout<<AcA[j]<<"\n";
            }
           // cout<<"BcA[j]:\n";
            for(size_t j=0;j<BcA.size();++j) {mol[NRA+BcA[j]] = oldMol; //cout<<BcA[j]<<"\n";
            } // A part retains oldMol number
            //cout<<"AcB[j]:\n";

            for(size_t j=0;j<AcB.size();++j) {mol[AcB[j]] = newMol;// cout<<"AcB[j]: "<<AcB[j]<<"mol["<<AcB[j]<<"]: "<<mol[AcB[j]]<<"\n";
            }
            //cout<<"BcB[j]:\n";
            for(size_t j=0;j<BcB.size();++j) {mol[NRA+BcB[j]] = newMol;// cout<<BcB[j]<<"\n";
            }

        }

    }
       // if(JA==328 && mA==606) {cout<<"reverse";system("pause");}
//cout<<"updated: mA"<<mol[JA]<<"mB"<<mol[NRA+JB]<<"\n";
}
int getJunctionDistAA(size_t J1,size_t J2)
{ // return -1 if J1 J2 are not connected, otherwise return minimal distance between them
// find minimum distance between J1 and neighbour of J2
    if(J1 == J2) return 0;
    if(neighA[J2].size()<1) return -1;
    if(JunctionDist[J1][neighA[J2][0]] == 0) return -1;
    size_t minD = MAX,JB; // Equivalent to: int minD = MAX; int JB;
    for(size_t i=0;i<neighA[J2].size();++i) {
            JB = neighA[J2][i];
            if(JunctionDist[J1][JB] < minD)
                minD = JunctionDist[J1][JB];
    }
    return minD+1;
}
int getJunctionDistBB(size_t J1,size_t J2)
{
    if(J1 == J2) return 0;
    if(neighB[J2].size()<1) return -1;
    if(JunctionDist[neighB[J2][0]][J1] == 0) return -1;
    size_t minD = MAX,JA;
    for(size_t i=0;i<neighB[J2].size();++i) {
        JA = neighB[J2][i];
        if(JunctionDist[JA][J1] < minD)
            minD = JunctionDist[JA][J1];
    }
    return minD+1;
}
size_t dist(unsigned char JunctionDist)
{
// JunctionDist should always be odd, since only different types of junctions react
    if(JunctionDist == 0) return 0;
    else return (JunctionDist+1)/(2*(dA+dB));
}
//Initializes random number generator with seed
//RNG is Mersenne Twister MT19937 algorithm
void RNG_initialize(unsigned long seed)
{
    mt[0]= seed & 0xffffffffUL;
    for(mti=1; mti<624; mti++){
    mt[mti] = (1812433253UL*(mt[mti-1]^(mt[mti-1] >> 30)) + mti);
    mt[mti] &= 0xffffffffUL;
    /* for >32 bit machines */
    }
    double dn = 3.442619855899;
    int i;
    const double m1 = 2147483648.0;
    double q;
    double tn = 3.442619855899;
    const double vn = 9.91256303526217E-03;
    q = vn/exp(-0.5*dn*dn);
    kn[0] = (dn/q)*m1;
    kn[1] = 0;
    wn[0] = q/m1;
    wn[127] = dn/m1;
    fn[0] = 1.0;
    fn[127] = exp(-0.5*dn*dn);
    for(i = 126; i >= 1; i--){
    dn = sqrt(-2*log(vn/dn + exp(-0.5*dn*dn)));
    kn[i+1] = (dn/tn)*m1;
    tn = dn;
    fn[i] = exp(-0.5*dn*dn);
    wn[i] = dn/m1;
    }
}
//Returns a random long between 0 and 4294967295
unsigned long rand32()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A for x=0,1 */
    if(mti >= 624){
    int kk;
    for(kk=0;kk<227;kk++){
    y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
    mt[kk] = mt[kk+397] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    for(;kk<623;kk++){
    y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
    mt[kk] = mt[kk-227] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    y = (mt[623]&UPPER_MASK)|(mt[0]&LOWER_MASK);
    mt[623] = mt[396] ^ (y >> 1) ^ mag01[y & 0x1UL];
    mti = 0;
    }
    y = mt[mti++];
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    //y=0.5*4294967296.0;
    return y;
}
