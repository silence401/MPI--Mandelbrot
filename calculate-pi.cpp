#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <time.h>
using namespace std;
#define  actual_pi  3.141592653589793238462643
int main(int argc,char** argv)
{
        int n;
        double s,e;
        double init_size,part_sum=0,x1,x2,temp_pi=0,cal_pi=0;
        int rank, flag=0, namelen;  
        int proc_num;
        char processor_name[MPI_MAX_PROCESSOR_NAME];   
        MPI_Status status;   
        MPI_Init(&argc,&argv);   
        MPI_Comm_size(MPI_COMM_WORLD,&proc_num);   
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);   
        MPI_Get_processor_name(processor_name,&namelen);   
        if(rank==0)
        {
            cout<<"输入切割的块数:";
            cin>>n;
            s=MPI_Wtime();
         //   MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
        }
         MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        init_size=1.0/n;
        for(int i=rank+1;i<=n;i+=proc_num)
        {
            x1=init_size*(i);
            x2=init_size*(i-1);
            x1=4/(1+x1*x1);
            x2=4/(1+x2*x2);
            part_sum=part_sum+x1+x2;
        }
        temp_pi=init_size*part_sum/2;
        MPI_Reduce(&temp_pi,&cal_pi,1,MPI::DOUBLE,MPI::SUM,0,MPI_COMM_WORLD);
        if(rank==0)
        {
              
        ///cout.flags(ios::fixed);  
        // cout.precision(4); //设置输出精度，  
      
            e=MPI_Wtime();
             cout<<"time: "<<e-s<<endl;
            cout<<"PI : "<<fabs(cal_pi)<<endl;
            cout<<"ERROR: "<<fabs(cal_pi - actual_pi)<<endl;

        }
     //else 
     //   {
    //        cout<<n<<endl;
//
     //   }
        MPI_Finalize();
}