/*

   fft_rla_for.c - run Fast Fourier Transform on RLA image
        steps :
            (1) read WF format image file composed of gaussian random B&W values (0-255)
            (2) perform forward FFT on image
            (3) write out magnitude file
            (4) write out phase file


    Jim Polk
    written: approx 1991

 */

#include <stdio.h>
#include <color_3D.h>
#include <fmt.h>
#include <imf.h>
#include <math.h>

typedef struct complex {
    float r,i;
} complex ;
typedef struct regular {
    float i;
} regular ;

#define PI 3.1415926535
#define NEWLINE 012
#define SPACE   040
#define PMODE   0755
#define READF   0
#define WRITF   1
#define ERR     -1


complex *compalloc(nn)
int nn;
{
    return ((complex *)malloc(nn*sizeof(complex)));
}

regular *regalloc(nn)
int nn;
{
    return ((regular *)malloc(nn*sizeof(regular)));
}




main (argc, argv)
int argc;
char *argv[];
{

    /* normal C stuff */

    complex *in_struct,*out_struct;
    regular *mag_struct,*phase_struct;
    complex *in_start, *out_start;
    regular *mag_start, *phase_start;
    register complex *sptr,*dptr;
    int isign;
    int n;
    float max,min,temp;
    float tempa,tempb,tempc;

    /* WF image read stuff */

    IMF_OBJECT          *imf_r_obj;
    int                 scan,scan2;                   /* current scanline */
    int                 pix,pix2;                    /* current pixel */
    POINTER             *line_buff;             /* pointer to line buffers */
    U_CHAR              *red, *green, *blue;    /* cur pixel pointer */

    /*   output stuff */

    int       magout,phaseout;    /* pointer to binary output files */
    int     siz;
    int        numbytes;

    /***************************/
    /*   open files            */
    /***************************/

    if (argc != 4)
    {
        fprintf(stderr,"usage:  rla_fft_for imagein.rla mag_out.mag phase_out.pha \n");
        exit (-1);
    }

    if ((imf_r_obj = IMF_open_read(argv[1])) == NULL)
    {
        fprintf (stderr, "imf_read: cannot open input image %s \n", argv[1]);
        exit (-1);
    }


    /********************************************/
    /* open the image file and write the header */
    /********************************************/


    if(( magout = creat(argv[2],PMODE)) == ERR )
    {
        fprintf(stderr,"can't open %s for output \n",argv[2]);
        exit(1);
    }


    if(( phaseout = creat(argv[3],PMODE)) == ERR )
    {
        fprintf(stderr,"can't open %s for output \n",argv[3]);
        exit(1);
    }




    /***********************************************/
    /*   allocate some memory to hold raster files */
    /***********************************************/

    n = 512;        /* image square size */
    isign = 1;        /* forward FFT */


    if(( in_struct  = compalloc( (n*2)*n) ) == NULL )
    {
        fprintf(stderr,"can't allocate in_struct \n");
        exit(1);
    }
    if(( out_struct = compalloc( (n*2)*n) ) == NULL )
    {
        fprintf(stderr,"can't allocate out_struct \n");
        exit(1);
    }
    if(( mag_struct = regalloc(n*n) ) == NULL )
    {
        fprintf(stderr,"can't allocate mag_struct \n");
        exit(1);
    }

    siz = (n*n) * sizeof(mag_struct);
    printf("size of mag_struct = %d \n",siz);

    if(( phase_struct = regalloc(n*n) ) == NULL )
    {
        fprintf(stderr,"can't allocate phase_struct \n");
        exit(1);
    }


    in_start = in_struct;
    out_start = out_struct;
    mag_start = mag_struct;
    phase_start = phase_struct;

    printf("\n     ...allocated memory \n");


    /*********************/
    /*  read input file  */
    /*********************/

    for (scan=0; scan<512; scan++)
    {

        if ((*imf_r_obj->scan)(imf_r_obj->data,scan,&line_buff) != IMF_C_NORMAL)
        {
            fprintf(stderr, "imf_read: read error on scan %d\n", scan);
        }

        red   = (U_CHAR *)line_buff[0];    /* input file */
        green = (U_CHAR *)line_buff[1];
        blue  = (U_CHAR *)line_buff[2];

        for (pix=0; pix<512; pix++)
        {
            tempa = -1;
            tempb = (scan+pix);

            tempc = pow(tempa,tempb);

            in_struct->r = *red++;
            in_struct->r = in_struct->r * tempc;
            in_struct->i = 0.0;
            in_struct++;
        }

    }                              /*  end for-each scan */

    (void)(*imf_r_obj->close)(imf_r_obj);         /* close input file */

    /****************************************************************************/
    /*   entire input image should now be contain in in_struct memory structure */
    /****************************************************************************/


    /* get pointer back to start */

    in_struct = in_start;

    /* transform each scan line */

    printf("\n     first fft...");

    for(scan=0; scan<512; scan++)
    {
        fft(n,isign,in_struct+(scan*n));
    }


    /* transpose rows and columns */

    transpose(n,in_struct,out_struct);


    /* transform each scan line */

    printf("\n     second fft...");

    for(scan=0; scan<512; scan++)
    {
        fft(n,isign,out_struct+(scan*n));
    }

    out_struct = out_start;


    /**************************************************/
    /* search complex number array for min/max values */
    /**************************************************/


    printf("\n     searching for min and max values ....");

    max = 0.0 ;
    min = 0.0 ;

    for(scan=0; scan<512; scan++)
    {
        for(pix=0; pix<512; pix++)
        {
            temp  = sqrt(( fabs(out_struct->r) * fabs(out_struct->r) ) + ( fabs(out_struct->i) * fabs(out_struct->i) ));

            if(temp > max ) max = temp;
            if(temp < min ) min = temp;

            out_struct++;
        }
    }

    printf("\n     max = %f   min = %f   ",max,min);

    out_struct = out_start;
    mag_struct = mag_start;
    phase_struct = phase_start;

    /******************************************/
    /* create magnitude and phase image files */
    /******************************************/

    printf("\n     calculating magnitude and phase...");

    for (scan=0; scan<512; scan++)
    {

        for (pix=0; pix<512; pix++)
        {
            mag_struct->i = sqrt(( fabs(out_struct->r) * fabs(out_struct->r) ) + ( fabs(out_struct->i) * fabs(out_struct->i) ));

            if(out_struct->r == 0)
                phase_struct->i = 0.0;
            else
                phase_struct->i = (float)(atan)( out_struct->i / out_struct->r);

            out_struct++;
            mag_struct++;
            phase_struct++;
        }

    }                /* end of for-each scan */

    in_struct = in_start;
    out_struct = out_start;
    mag_struct = mag_start;
    phase_struct = phase_start;


    printf("\n     writing magnitude & phase files...");

    /* magnitude image */

    siz = (n*n) * sizeof(mag_struct);

    printf("\n     siz of mag_struct = %d \n",siz);

    if(( numbytes = write(magout,mag_struct,siz)) != siz)
        fprintf(stderr,"Error writing out magnitude structure \n");

    /* phase image */

    siz = (n*n) * sizeof(phase_struct);

    printf("\n     siz of phase_struct = %d \n",siz);

    if(( numbytes = write(phaseout,phase_struct,siz)) != siz)
        fprintf(stderr,"Error writing out phase structure \n");

    printf("done\n");

    /*****************************/
    /*     close the file        */
    /*****************************/

    free(in_start);
    free(out_start);
    free(mag_start);
    free(phase_start);

    close(magout);
    close(phaseout);


}                /********* end of main ********/



fft(n,isign,data)
int n,isign;
register complex data[];
{

    register int i;
    int j, k, ln, nv2;
    complex temp, u, w;
    int le, le1, l;

    powercheck(n);

    if(isign > 0)
        isign = 1;
    else
        isign = -1;

    for(ln=0, i=n; i>1; i>>=1)
        ln++;

    nv2 = n / 2;
    j=0;
    n--;
    for(i=0; i<n; i++)
    {
        if(i < j)
        {
            temp = data[j];
            data[j] = data[i];
            data[i] = temp;
        }
        k = nv2;
        while(j >= k)
        {
            j -= k;
            k >>= 1;
        }
        j += k;
    }

    n++;
    le = 2;
    le1 = 1;
    for(l=0; l<ln; l++)
    {
        register float tempr, tempi;
        register float ur,ui;
        register float wr,wi;
        register complex *ipptr, *iptr;

        ur = 1.0;
        ui = 0.0;
        wr =  cos(isign * PI / le1);
        wi = -sin(isign * PI / le1);
        for(j=0; j<le1; j++)
        {
            for(i=j; i<n; i+=le)
            {
                iptr = data + i;
                ipptr = iptr + le1;
                tempr = ur * ipptr->r - ui * ipptr->i;
                tempi = ur * ipptr->i + ui * ipptr->r;
                ipptr->r = iptr->r - tempr;
                ipptr->i = iptr->i - tempi;
                iptr->r += tempr;
                iptr->i += tempi;
            }
            tempr = ur;
            ur = wr*ur - wi*ui;
            ui = wr*ui + wi*tempr;
        }
        le <<= 1;
        le1 <<= 1;
    }
    if(isign == 1);
    {
        for(i=0; i<n; i++)
        {
            data[i].r /= n;
            data[i].i /= n;
        }
    }
}


powercheck(n)
int n;
{
    int i;

    for(i=0; i<32; i++)
    {
        if(n & 1)
            break;
        n >>= 1;
    }
    if( n != 1)
    {
        fprintf(stderr,"fft: n points must be a power of 2. \n");
        exit(1);
    }
}


transpose(n,in_struct,out_struct)
register int n;
complex *in_struct;
complex *out_struct;
{

    register int scan,pix;
    register complex *sptr, *dptr;

    printf("\n     transposing image...");
    sptr = in_struct;

    for(scan=0; scan<512; scan++)
    {
        dptr = out_struct + scan;

        for(pix=0; pix<512; pix++)
        {
            *dptr = *sptr++;
            dptr += n;
        }

    }


}






