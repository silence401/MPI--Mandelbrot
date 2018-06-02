    #include <mpi.h>   
    #include <stdio.h>   
    #include <GL/glut.h>   
    #include <time.h>   
    #include <stdlib.h>   
    #include <unistd.h>
    #include <stdlib.h>
    //定义status返回的MPI_TAG   
    #define Over 0   
    #define Redo 1   
    #define Exit 2   
    const double MIN_INTERVAL = 1e-10;   
    const double MAX_MAGNITUDE = 4.0;   
    const short  MAX_ITERATE_DEPTH = 20000;   
    const short  MAX_ZOOM_DEPTH = 500;   
    const int    MAX_WIN_SIZE = 800;   
    const int    PIXEL_NUM = 800*800;   
    //定义复数结构体   
    struct ComplexNumber   
    {   
        double real;   
        double imag;   
    };   
    //定义绘制区域结构体   xloop,yloop 表示整个窗口区域要绘制的像素值sub_xloop,sub_yloop表示每一个进程
    struct Area   
    {   
        double  Xmin;   
        double  Ymin;   
        double  Xarea;   
        double  Yarea;   
        double  dx;   
        double  dy;   
        int     xloop;   
        int     yloop;   
        int     sub_xloop;   
        int     sub_yloop;   
    };   
    int win_size;  
    int count=0;
    int zoom_depth;   
    Area infos[MAX_ZOOM_DEPTH];   
    Area part[10000];
    Area info;   
    short recv[MAX_WIN_SIZE][MAX_WIN_SIZE];   
    short indices[MAX_WIN_SIZE][MAX_WIN_SIZE];   
    double color[MAX_ITERATE_DEPTH+1][3];   
    int proc_num, slave_num;   
    int master_exit;    
    MPI_Datatype AreaType;    
    double start_time, end_time;    
    void Redraw(int w, int h);   
    void Mouse(int button, int state, int x, int y);   
    void Display();   
    void Idle();   
    void Unproject(int x, int y, double *objx, double *objy);
    void ResetWin();  
    //迭代函数   
    ComplexNumber f(ComplexNumber z, ComplexNumber c)   
    {   
        ComplexNumber result;   
        result.real=c.real+z.real*z.real-z.imag*z.imag;   
        result.imag=c.imag+z.imag*z.real+z.real*z.imag;   
        return result;  
    }   
    void Mandelbrot(double Xmin,double dx, int xloop, int xfrom, double Ymin,  double dy, int yloop, int yfrom)   
    {   
        int x, y, k;   
        ComplexNumber c, z;   
        for (x=0; x<xloop; x++)   
        {   
            c.real = Xmin+x*dx;   
            for (y=0; y<yloop; y++)   
            {   
                c.imag = Ymin+y*dy;   
                z.real = z.imag = 0.0f;   
                k = 0;   
                while (k<MAX_ITERATE_DEPTH && (z.real*z.real+z.imag*z.imag)<=MAX_MAGNITUDE)   
                {   
                    z = f(z, c);   
                    k++;   
                }   
                recv[xfrom+x][yfrom+y] = k;   
            }   
        }   
    }   
    void NewArea(double x1, double y1, double x2, double y2)   
    {   
        double t;   
        if (x1 > x2)    
        {   
            t = x1;   
            x1 = x2;   
            x2 = t;   
        }   
        if (y1 > y2)    
        {   
            t = y1;   
            y1 = y2;   
            y2 = t;   
        }   
        if (zoom_depth > 0)   
        {   
            info = infos[zoom_depth-1];   
            if (x1 < info.Xmin)   
                x1 = info.Xmin;      
            if (x2 > info.Xmin+info.Xarea)   
                x2 = info.Xmin+info.Xarea;   
            if (y1 < info.Ymin)   
                y1 = info.Ymin;   
            if (y2 > info.Ymin+info.Yarea)   
                y2 = info.Ymin+info.Yarea;   
        }   
        info.Xmin = x1;   
        info.Ymin = y1;   
        info.Xarea = x2-x1;   
        info.Yarea = y2-y1;   
        if (info.Xarea>=info.Yarea)   
        {   
            info.xloop = MAX_WIN_SIZE;   
            info.yloop = (int)(info.Yarea/info.Xarea*MAX_WIN_SIZE);   
        }   
        else   
        {   
            info.yloop = MAX_WIN_SIZE;   
            info.xloop = (int)(info.Xarea/info.Yarea*MAX_WIN_SIZE);   
        }      
        info.sub_xloop = info.xloop/slave_num;   
        info.sub_yloop = info.yloop/slave_num;
        info.xloop = info.sub_xloop*slave_num;   
        info.yloop = info.sub_yloop*slave_num;   
        info.dx = info.Xarea/info.xloop;   
        info.dy = info.Yarea/info.yloop;   
        infos[zoom_depth++] = info;   
        for(int i=0;i<slave_num*slave_num;i++)
        {
            part[i].xloop=i;
            part[i].sub_xloop=info.sub_xloop;
            part[i].sub_yloop=info.sub_yloop;
            part[i].Xarea=info.Xarea;
            part[i].Yarea=info.Yarea;
            part[i].Xmin=info.Xmin;
            part[i].Ymin=info.Ymin;
            part[i].dx=info.dx;
            part[i].dy=info.dy;
        }
        count=slave_num-1;
    }   
    //向每个从机发送新的绘制区域信息，MPI_TAG=Redo   
    void BroadcastRedo()   
    {   
       for (int i=1; i<proc_num; i++)   
        MPI_Send(&part[i-1],1,AreaType,i,Redo,MPI_COMM_WORLD);   
        srand((unsigned int)time(NULL));
        start_time = MPI_Wtime();   
        ResetWin();   
    }   
    //向每个从机发送退出消息，MPI_TAG=Exit   
    void BroadcastExit()    
    {   
        master_exit = 1;   
        for (int i=1; i<proc_num; i++)   
            MPI_Send(&master_exit,1,MPI_INT,i,Exit,MPI_COMM_WORLD);   
        MPI_Finalize();   
    }   
    //初始化OpenGL窗体信息   
    void ShowWin()   
    {       
        int i; 
        // glutInit(&argc,argv);  
        glutInitDisplayMode(GLUT_DOUBLE);   
        glutInitWindowSize(640,640);   
        glutInitWindowPosition(100,100);   
        glutCreateWindow("Mandelbrot");   
        srand((unsigned int)time(NULL));   
        for (i=0; i<MAX_ITERATE_DEPTH; i++)   
        {   
            color[i][0] = (GLdouble)rand()/(GLdouble)RAND_MAX;   
            color[i][1] = (GLdouble)rand()/(GLdouble)RAND_MAX;   
            color[i][2] = (GLdouble)rand()/(GLdouble)RAND_MAX;   
        }      
        NewArea(-2.0,-2.0,2.0,2.0);    
        BroadcastRedo();    
        glClearColor(0.0, 0.0, 0.0, 1.0);   
        glutDisplayFunc(Display);   
        glutReshapeFunc(Redraw);   
        glutMouseFunc(Mouse);   
        glutIdleFunc(Idle);   
        glutMainLoop();  
    }   
    // OpenGL窗口显示函数   
    void Display()    
    {   
        //printf("1");
        int i, j;   
        glClear(GL_COLOR_BUFFER_BIT);     
        glBegin(GL_POINTS);   
        for (i=0; i<info.xloop; i++)   
        {   
            for (j=0; j<info.yloop; j++)   
            {   
                if (MAX_ITERATE_DEPTH != indices[i][j])    
                {   
                    glColor3dv(color[indices[i][j]]);   
                    glVertex2d(info.Xmin+i*info.dx, info.Ymin+j*info.dy);      
                }       
              
            }   
           
            usleep(500);
        }   
        glEnd();     
        glFlush();   
        glutSwapBuffers();  
       // printf("2");  
    }   
    // OpenGL窗口调整事件响应函数   
    void Redraw(int w, int h)    
    {   
        if (w <= h)   
        {   
            win_size = w;   
            glViewport(0, (h-w)/2, w,w);   
        }   
        else   
        {   
            win_size = h;   
            glViewport((w-h)/2, 0, h,h);   
        }   
    }    
    //OpenGL窗口闲置函数   
    void Idle()    
    {   
        //printf("34");
        static int recved = 0;   
        int pos;
        int i, j, x, y;   
        int xfrom, yfrom;   
        int flag, slave_rank;   
        MPI_Status status;   
        MPI_Iprobe(MPI_ANY_SOURCE,Over,MPI_COMM_WORLD, &flag, &status);   
        if (!flag)  return;   
        MPI_Recv(&pos,1,MPI_INT,status.MPI_SOURCE,Over,MPI_COMM_WORLD,&status);
        MPI_Recv(recv,PIXEL_NUM,MPI_SHORT,status.MPI_SOURCE,Over,MPI_COMM_WORLD,&status);   

        if (slave_num*slave_num== ++recved)   
        {   
            recved = 0;   
            end_time = MPI_Wtime();    
            printf("wall clock time = %f\n", end_time-start_time);   
        }   
        slave_rank = status.MPI_SOURCE-1;   
        printf("recieve from slave %d\n", slave_rank);   
        i = slave_rank;   
        xfrom = (pos%slave_num)*info.sub_xloop;   
        yfrom = (pos/slave_num)*info.sub_yloop;       
        for(x=0; x<info.sub_xloop; x++)   
        {   
                for (y=0; y<info.sub_yloop; y++)   
                    indices[xfrom+x][yfrom+y] = recv[xfrom+x][yfrom+y];   
        }
        if(count<slave_num*slave_num)
       {
         //  printf("%d\n",count);
       //  printf("here\n");
        MPI_Send(&part[count],1,AreaType,slave_rank+1,Redo,MPI_COMM_WORLD);   
        count++;
        }
        glutPostRedisplay();
    }   
    //鼠标监听
void Mouse(int button, int state, int x, int y)   
    {   
        static int x1, y1;   
        if (GLUT_LEFT_BUTTON == button)   
        {   
            if (GLUT_DOWN == state)  
            {   
                x1 = x;   
                y1 = y;   
            }   
            if (GLUT_UP == state)   
            {   
              //  if (zoom_depth == MAX_ZOOM_DEPTH)   
                 //   return;   
                double objx1,objy1,objx2,objy2;   
                Unproject(x1, y1, &objx1, &objy1);    
                Unproject(x, y, &objx2, &objy2);   
                NewArea(objx1,objy1,objx2,objy2);   
                BroadcastRedo();   
                glutPostRedisplay();   
            }   
        } 
    } 
            // 把窗口坐标(X,Y)转换成对象坐标   
    void Unproject(int x, int y, double *objx, double *objy)    
    {   
        GLdouble modelMatrix[16], projMatrix[16];   
        GLint viewport[4];   
        GLfloat winz;   
        GLdouble objz;   
        glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);   
        glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);   
        glGetIntegerv(GL_VIEWPORT, viewport);      
        glReadPixels(x,y,1,1,GL_DEPTH_COMPONENT,GL_FLOAT, &winz);   
        gluUnProject((GLdouble)x,(GLdouble)(win_size-y),(GLdouble)winz,modelMatrix,projMatrix,viewport,objx,objy,&objz);   
    }   
    void ResetWin()   
    {   
        glMatrixMode(GL_PROJECTION);   
        glLoadIdentity();   
        info = infos[zoom_depth-1];   
        if (info.Xarea >= info.Yarea)   
        {   
           printf("%lf %lf %lf %lf\n",info.Xmin,info.Ymin,info.Xmin+info.Xarea,info.Ymin+info.Yarea);
           gluOrtho2D(info.Xmin,info.Xmin+info.Xarea,info.Ymin-(info.Xarea-info.Yarea)/2,info.Ymin+info.Yarea+(info.Xarea-info.Yarea)/2);   
           //   gluOrtho2D(info.Xmin,info.Xmin+info.Xarea,info.Ymin+info.Xarea,info.Ymin);
            // printf("%lf %lf %lf %lf\n",info.Xmin,info.Xmin+info.Xarea,info.Ymin-(info.Xarea-info.Yarea)/2,info.Ymin+info.Yarea+(info.Xarea-info.Yarea)/2);
        }   
        else   
        { 
            printf("%lf %lf %lf %lf\n",info.Xmin,info.Ymin,info.Xmin+info.Xarea,info.Ymin+info.Yarea);
            gluOrtho2D(info.Xmin-(info.Yarea-info.Xarea)/2,info.Xmin+info.Xarea+(info.Yarea-info.Xarea)/2,info.Ymin,info.Ymin+info.Yarea); 
           //  gluOrtho2D(info.Xmin,info.Xmin+info.Xarea,info.Ymin+info.Xarea,info.Ymin);
             //printf("%lf %lf %lf %lf\n",info.Xmin,info.Xmin+info.Xarea,info.Ymin-(info.Xarea-info.Yarea)/2,info.Ymin+info.Yarea+(info.Xarea-info.Yarea)/2);
          //  gluOrtho2D(info.Xmin,info.Ymin,info.Xarea,info.Yarea);  
        }   
        glMatrixMode(GL_MODELVIEW);   
        glLoadIdentity();   
        int i, j, xloop, yloop;   
        xloop = infos[zoom_depth-1].xloop;   
        yloop = infos[zoom_depth-1].yloop;   
        for (i=0; i<xloop; i++)   
        {   
            for (j=0; j<yloop; j++)   
                indices[i][j] = MAX_ITERATE_DEPTH;   
        }   
        printf("//***************************************//\n");   
        printf("//current depth = %d\n", zoom_depth);   
        printf("//Xmin=%f Xarea=%f Ymin=%f Yarea=%f\n",info.Xmin,info.Xarea,info.Ymin, info.Yarea);   
        printf("//dx=%f dy=%f\n", info.dx,info.dy);   
        printf("//xloop=%d yloop=%d subx=%d suby=%d\n",info.xloop, info.yloop,info.sub_xloop, info.sub_yloop);   
        printf("//***************************************//\n");   
    }   
    int main(int argc, char * argv[])   
    {   
        int rank, flag=0, namelen;   
        char processor_name[MPI_MAX_PROCESSOR_NAME];   
        MPI_Status status;   
        MPI_Init(&argc,&argv);   
        MPI_Comm_size(MPI_COMM_WORLD,&proc_num);   
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);   
        MPI_Get_processor_name(processor_name,&namelen);   
        MPI_Datatype types[2] = {MPI_DOUBLE,MPI_INT};   
        int blocklens[2]= {6,4};   
        MPI_Aint disp[2];   
        MPI_Address(&info, disp);   
        MPI_Address(&(info.xloop), disp+1);   
        disp[1] -= disp[0];   
        disp[0] = 0;   
        MPI_Type_struct(2, blocklens, disp, types, &AreaType);   
        MPI_Type_commit(&AreaType);   
        if (0 == rank)   
            fprintf(stderr,"Master Process  on %s\n", processor_name);   
        else   
            fprintf(stderr,"Slave Process %d on %s\n", rank-1, processor_name);   
        fflush(stderr);   
        MPI_Barrier(MPI_COMM_WORLD);   
       
        slave_num = proc_num-1;   
        if (0 == slave_num)   
        {   
            printf("WARNING:This program need at least one slave processor!\n");   
            return 0;   
        }   
        master_exit = 0;   
           
        if (0 == rank)   
        { 
            glutInit(&argc,argv);
            atexit(BroadcastExit);   
            ShowWin();   
        }   
        else   
        {   
           while (!flag)   
           {   
                MPI_Iprobe(0,MPI_ANY_TAG,MPI_COMM_WORLD, &flag, &status);   
                if (flag)   
                {   
                     //MPI_Recv(&info,1,AreaType,0,Redo,MPI_COMM_WORLD,&status);
                    if (Exit == status.MPI_TAG)   
                    {   
                        MPI_Finalize();   
                        return 0;   
                    }   
                if (Redo == status.MPI_TAG)   
                    {   
                        MPI_Recv(&info,1,AreaType,0,Redo,MPI_COMM_WORLD,&status);   
                        printf("%d\n",info.xloop);
                         //Mandelbrot(double Xmin,double dx, int xloop, int xfrom, double Ymin,  double dy, int yloop, int yfrom)  
                        double txmin=info.Xmin+(info.xloop%slave_num)*info.Xarea/slave_num;
                        int txfrom=(info.xloop%slave_num)*info.sub_xloop;
                        double tymin=(info.xloop/slave_num)*info.Yarea/slave_num+info.Ymin;
                        int tyfrom=(info.xloop/slave_num)*info.sub_yloop;
                      //  Mandelbrot(info.Xmin+(info.xloop%slave_num)*info.Xarea/slave_num, info.dx,info.sub_xloop, (info.xloop*info.sub_xloop,info.Ymin+i*info.Yarea/slave_num, info.dy,info.sub_yloop, i*info.sub_yloop);   
                        Mandelbrot(txmin,info.dx,info.sub_xloop,txfrom,tymin,info.dy,info.sub_yloop,tyfrom);
                        MPI_Send(&info.xloop,1,MPI_INT,0,Over,MPI_COMM_WORLD);
                        MPI_Send(recv,PIXEL_NUM,MPI_SHORT,0,Over,MPI_COMM_WORLD);   
                        printf("send from slave %d\n", rank-1);   
                    }   
                //}   
                } 
                flag=0;    
           }
         //  while(1);
        }      
        return 0;   
    }  
