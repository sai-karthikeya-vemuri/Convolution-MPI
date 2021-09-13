#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <mpi.h>

class Convolution
{
private:
std::vector<std::vector<int>> image;
std::vector<std::vector<float>> kernel;
std::vector<std::vector<int>> output;

int rows;
int columns;

public:

Convolution(std::vector<std::vector<int>> image_input,std::vector<std::vector<float>> kernel_input,int rows,int columns)
  : rows(rows),columns(columns)
{


image = image_input;
kernel= kernel_input;
output.resize(rows, std::vector<int>(columns,0));
}
int index_number(int a, int max)
{
    int index;
    if (a < 0){
        index = max + a;
    }
    else if (a >= max){
        index =  a-max;
    }
    else{
        index = a;
    }
    return index;
}
std::vector<std::vector<int>> Convolute_Serial()
{
int indexr,indexc,a,b;
float sum;






//std::cout<<"Go ahead, program the convolution";
for(int x = 0; x<rows;++x)
 {
   for(int y=0; y<columns;++y) 
   {
     sum=0;
     for(int i=0;i<3;++i)
     {
     for(int j=0;j<3;j++)
      {
      a = x +i -1;
      b = y+j -1;
      indexr= index_number(a,rows);
      indexc= index_number(b,columns);
      sum += image[indexr][indexc] * kernel[i][j];
      
      
      }
     }
   if(sum>255) {sum=255;} else if(sum<0){sum=0;} 
   output[x][y]=sum;
   }
   
 } 
 return output;
}

std::vector<std::vector<int>> Convolute_Parallel(int rank,int p)
{
    //int rank;
    //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //int p;
    //MPI_Comm_size(MPI_COMM_WORLD, &p);
	// MPI_Barrier(MPI_COMM_WORLD);
	
	if(p==1 || p>rows){
	        std::cout<<"Serial done";
		return Convolute_Serial();}
        else
	{
	std::cout<<"In else \n";
	int rows_p = rows/p;
	int temp_array[rows][columns];
	for (int i =0; i < rows; ++i)
    {
        for (int j = 0; j < columns; ++j) {
            temp_array[i][j] = image[i][j];
        }
    }
	std::cout<<"Temp array copied succesfully \n";
	int dest1, dest2, src1, src2;                 
    int recHalo1[columns], recHalo2[columns];      
    int sendHalo1[columns], sendHalo2[columns];       
    std::vector<std::vector<int>> matrixHalo;  
	int matrixHaloArr[rows_p][columns];     
    float convProcMatrix[rows_p][columns];   
    float FinalMatrix [rows][columns]; 
	
    MPI_Scatter(temp_array, rows_p*columns, MPI_INT, 
        &matrixHaloArr, rows_p*columns, MPI_INT, 0, MPI_COMM_WORLD);
        
	std::cout<<"Done scattering \n";
	MPI_Request requestHandle[2];
        
       
    for (int i = 0; i < columns; i++){
        sendHalo1[i] = matrixHaloArr[0][i];
        sendHalo2[i] = matrixHaloArr[rows_p - 1][i];
    }    
    std::cout<<"Done Halo\n";
    if (p == 1) {//if there is only 1 process
        dest1 = 0;
        dest2 = 0;
        src1 = 0;
        src2 = 0;
    } 
    else {//if more than one process are there
        if (rank == 0) {
            dest1 = p - 1;
            dest2 = rank + 1;
            src1 = p - 1;
            src2 = rank + 1;
            } else if (rank == p - 1) {
                dest1 = rank - 1;
                dest2 = 0;
                src1 = rank - 1;
                src2 = 0;
            } else {
                dest1 = rank - 1;
                dest2 = rank + 1;
                src1 = rank - 1;
                src2 = rank + 1;
            }
	
    
	
	}
    MPI_Isend(&sendHalo1, columns, MPI_INT, dest1, 111, MPI_COMM_WORLD, &requestHandle[0]);
    MPI_Isend(&sendHalo2, columns, MPI_INT, dest2, 112, MPI_COMM_WORLD, &requestHandle[1]);
    MPI_Irecv(&recHalo1, columns, MPI_INT, src1, 112, MPI_COMM_WORLD, &requestHandle[1]);
    MPI_Irecv(&recHalo2, columns, MPI_INT, src2, 111, MPI_COMM_WORLD, &requestHandle[0]);
    MPI_Waitall(2, requestHandle, MPI_STATUSES_IGNORE);
        
	std::cout<<"Done isend and irecv\n";
	
	
	matrixHalo.resize(rows_p + 2, std::vector<int>(columns, 0));
    //     //convProcMatrix.resize(rows_p, std::vector<int>(columns, 0));
    for(int row=0;row<rows_p+2;row++){      
	    for (int i = 0; i < columns; i++){
            if(row==0){
	            matrixHalo[0][i] = recHalo1[i];
	        }   
	        else if(row==rows_p+1){
                matrixHalo[rows_p+1][i] = recHalo2[i];
	        }
	        else{
	            matrixHalo[row][i]=matrixHaloArr[row-1][i];
	        }
        }
	}
	int a,b,indexr,indexc;
	float sum;
	for(int x = 1; x<rows_p+1;++x)
    {
        for(int y=0; y<columns;++y) 
        {
        sum=0;
        for(int i=0;i<3;++i)
            {
                for(int j=0;j<3;++j)
                {
                    a = x +i -1;
                    b = y +j -1;
                    indexr= index_number(a,rows);
                    indexc= index_number(b,columns);
                    sum += matrixHalo[indexr][indexc] * kernel[i][j];
                }
            }
            if(sum>255) {sum=255;} else if(sum<0){sum=0;} convProcMatrix[x-1][y]=sum;
        }
    } 
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout<<"Done convolution process \n";
    MPI_Gather(&convProcMatrix, rows_p*columns, MPI_FLOAT, FinalMatrix,rows_p*columns , MPI_FLOAT, 0, MPI_COMM_WORLD);
    if(rank==0)
    {
        std::cout<<"Done gathering\n";
	    for(int i = 0;i<rows;++i){
	        for(int j=0;j<columns;++j){
	            output[i][j] = FinalMatrix[i][j];
	        }
	    } 
	}
    MPI_Barrier(MPI_COMM_WORLD);
    return output;
    } 
}

void create_ofile(std::string& ofilename)
{
    std::ofstream  outfile(ofilename);
    if(outfile.is_open())
    {
        outfile<<"P2"<<"\n";
        outfile<<rows<<" "<<columns<<"\n";
        outfile<<255<<"\n";
        for(int i =0;i<rows;i++){
            for (int j = 0;j<columns; j++){
                outfile<<output[i][j]<<" ";
            }
            outfile<<"\n";
        }
    }
    outfile.close();
}

~Convolution() {}

};

int main(int argc,char **argv)
{
  	MPI_Init(&argc,&argv);
	int rank;
    int pcom;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &pcom);
    int row = 0, col = 0, numrows = 0, numcols = 0, greyval = 0;
    std::ifstream infile("512.pgm");
    std::stringstream ss;
    std::string inputLine = "";
    // First line : version
    getline(infile,inputLine);
    if(inputLine.compare("P2") != 0) std::cerr << "Version error" << "\n";
    else std::cout << "Version : " << inputLine << "\n";


    // Continue with a stringstream
    ss << infile.rdbuf();
    // Third line : size
    ss >> numcols >> numrows >> greyval;
    std::cout << numcols << " columns and " << numrows << " rows" << "\n";

    std::vector<std::vector<int>> array;
    array.resize(numrows, std::vector<int>(numcols,0));
    // Following lines : data
    for(row = 0; row < numrows; ++row)
        for (col = 0; col < numcols; ++col) ss >> array[row][col];


    infile.close(); 
    std::vector<std::vector<float>> blur
    {
        {0.0625,0.125,0.0625},
        {0.125,0.25,0.125},
        {0.0625,0.125,0.0625}
    };
    std::vector<std::vector<float>> edge
    {
        {-1,-1,-1},
        {-1,8,-1},
        {-1,-1,-1}
    };
    std::vector<std::vector<float>> sharp
    {
        {0,-1,0},
        {-1,5,-1},
        {0,-1,0}
    };
    std::vector<std::vector<float>> identity
    {
        {0,0,0},
        {0,1,0},
        {0,0,0}
    };
    //Task 1:Edge

    Convolution C(array,edge,numrows,numcols);

    array=C.Convolute_Parallel(rank,pcom);
    if(rank == 0){
        std::string oname="512_parallel_edge.pgm";
        C.create_ofile(oname);
    }



/*
//Task2 : 5x Blur
for(int k =0;k<5;k++)
{
Convolution C(array,blur,numrows,numcols);

array=C.Convolute_Serial();
if(k==4)
{
Convolution C(array,sharp,numrows,numcols);
array=C.Convolute_Serial();
Convolution C1(array,edge,numrows,numcols);
array=C1.Convolute_Serial();
std::string oname="512pix_5xBlur_Edge_Sharp.pgm";
C1.create_ofile(oname);
}
}
*/
    MPI_Finalize();
    return 0;
}
