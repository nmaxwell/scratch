    
#include <mathlib/math/std_math.h>

#include <mathlib/math/transforms/fft.h>

#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>

//#include <pthread.h>




double get_time_resolution()
{
    double T0=get_real_time();
    for (int k=0; k<1000000; k++)
    {
        double t1= get_real_time()-T0;
        double t2 = get_real_time()-T0;
        
        if ( t2-t1 != 0)
            return (t2-t1);
    }
}









void *threadFunc(void *args)
{
    double T0 = *(double*)args;
    int tres = *(int*)((char*)args+8);    
    
    for (int k=0; k<100; k++)
    {
        cout << "t\t" << (get_real_time()-T0)*tres << endl;
    }
    
    
    
	return NULL;
}


int main()
{
    std_setup();
    
    

    double T0;
    int tres=0;
    tres = 1.0/get_time_resolution();
    T0=get_real_time();
    cout << tres << endl;
    
    
    
    
    pthread_t pth;
    
    char args[100];
    *(double*)args =T0;
    *(int*)(args+8)=tres;
    
    
    pthread_create( &pth, NULL, threadFunc, (void*)args );
    
    pthread_join(pth, NULL);
    
    for (int k=0; k<100; k++)
    {
        cout << "m\t" << (get_real_time()-T0)*tres << endl;
    }
    
    
    
	
    
    std_exit();
}










/*-
    int n= pow(2,12);
    double * u=ml_alloc<double> (n);
    complex<double> * v=ml_alloc<complex<double> >  (n);
    
    fft(u,v,n);
    */



