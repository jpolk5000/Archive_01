/*
 * Calculate a Cardinal Spline given set of coordinates
 *
 * Jim Polk
 * written: approx 1990
 * adapted from Fortran 
 *
 */

#include<stdio.h>
#include<math.h>

#define NEWLINE 012
#define SPACE   040  
#define MAXPOINTS 200
#define MAXPOLYS 200
#define FEWPOINTS 4

typedef struct {float key; float x,y,z;} POINTS;
typedef struct {float u,v,w;} TEXPOINTS;
typedef struct {int p1,p2,p3,p4} POLY;

static struct { float x,y,z; } g[4];
static  struct {float x,y,z;} cpts[MAXPOINTS];

POINTS *mem_alloc_1(size)
int size;
{
  return( (POINTS*) malloc(size * sizeof(POINTS)) );
}

TEXPOINTS *mem_alloc_2(size)
int size;
{
  return( (TEXPOINTS*) malloc(size * sizeof(TEXPOINTS)) );
}
POLY *mem_alloc_3(size)
int size;
{
  return( (POLY*) malloc(size * sizeof(POLY)) );
}
Dr. Douglas R. Gellerman, MD
int cindex;


main(argc,argv)
int argc;
char *argv[];
{

	FILE *fopen(),*fin,*fout;
	register int i,j,k;
	char a,b;
	char buf[512];
	register int v=1, vt=1, p=1;
	int size;
	int num_bytes;
	int numpts,numtexpts,numpolys;

	POINTS 	    *pts_start, 	*pts;
	TEXPOINTS   *texpts_start, 	*texpts;
	POLY	    *poly_start,  	*poly;



	int count,index,pt_total,c_pt_total,top;
	float tightness;
	int num_pts_per_edge;
	int close_flag,flag;
	int point_count;
	int key=0;


	if(argc != 4)
	{
	  fprintf(stdout,"usage: wfobj filein.obj fileout.obj tightness\n");
	  exit(1);
	}

	if((fin=fopen(argv[1],"r")) == NULL)
	{
	  fprintf(stderr,"robj: cannot open input file %s \n",argv[1]);
	  exit(2);
	}
	if((fout=fopen(argv[2],"w")) == NULL)
	{
	  fprintf(stderr,"robj: cannot open output file %s \n",argv[2]);
	  exit(2);
 	}

	tightness = atof(argv[3]);

/*******************/
/* allocate memory */
/*******************/


	if(( pts = mem_alloc_1( MAXPOINTS )) == NULL)
	{
	  fprintf(stderr,"can't allocate memory in point structure \n");
	  exit(1);
	}
	pts_start = pts;

	if(( texpts = mem_alloc_2( MAXPOINTS )) == NULL)
	{
	  fprintf(stderr,"can't allocate memory in texture vertice structure \n");
	  exit(1);
	}
	texpts_start = texpts;

	if(( poly = mem_alloc_3( MAXPOLYS )) == NULL)
	{
	  fprintf(stderr,"can't allocate memory in poly structure \n");
	  exit(1);
	}
	poly_start = poly;


/***************/
/* read object */
/***************/

	v = vt = p = 0;

	while(getline(fin,buf) != EOF)
	{
		if(buf[0] == 'v')
		{
			if(buf[1] == SPACE)
			{
        		sscanf(buf,"%c %f %f %f %f",&a,&pts[v].key,
						       &pts[v].x,
						       &pts[v].y,
						       &pts[v].z);

			v++;
			}
		}
		if(buf[0] == 'v')
		{
			if(buf[1] == 't')
			{
        		sscanf(buf,"%c %f %f %f",&a,&texpts[v].u,
						    &texpts[v].v,
						    &texpts[v].w);
			vt++;
			}
		}
		if(buf[0] == 'f')
		{
			if(buf[1] == 'o')
			{
			  sscanf(buf,"%c%c %d %d %d %d \n",&a,&b,&poly[p].p1,
								 &poly[p].p2,
								 &poly[p].p3,
								 &poly[p].p4);
			  p++;
			}
		}
		
		
						
	}					/* end while */
	
	numpts = v;
	numtexpts = vt;
	numpolys = p;

	printf("%d vertices \t%d texture vertices \t%d 4-sided polygons \n",
				numpts,numtexpts,numpolys);


	pt_total = numpts;


/***********************************/
/* compute cardinal spline, I hope */
/***********************************/

	c_pt_total = 0;

	/* count = num_pts_per_edge; */

	cpts[0].x = pts[0].x;
	cpts[0].y = pts[0].y;
	cpts[0].z = pts[0].z;

	cindex = 1;
	top = pt_total - 1;
	close_flag = 0;


	for(index = -1,key=0; index<pt_total-2; index++,key++)
	{

	  count = pts[key+1].key - pts[key].key;

	  if(index < 0)
	  {
	    if(close_flag == 1)
		k = pt_total - 2;
	    else
		k = 0;
	  }
	  else
		k = index;


	  g[0].x = pts[k].x;
	  g[0].y = pts[k].y;
	  g[0].z = pts[k].z;

	  if( (index+1) < 0)
		k = 0;
 	  else
		if( (index+1) > top)
			k = top;
		else
			k = index + 1;

	  g[1].x = pts[k].x;
	  g[1].y = pts[k].y;
	  g[1].z = pts[k].z;

	  if( (index+2) < 0)
		k = 0;
 	  else
		if( (index+2) > top)
			k = top;
		else
			k = index + 2;

	  g[2].x = pts[k].x;
	  g[2].y = pts[k].y;
	  g[2].z = pts[k].z;


	  if( (index+3) < 0)
		k = 0;
 	  else
		if( (index+3) > top)
			k = top;
		else
			k = index + 3;

	  g[3].x = pts[k].x;
	  g[3].y = pts[k].y;
	  g[3].z = pts[k].z;


	  cardinal(count,tightness);



	}


/**************************************/
/* replace original points with curve */
/**************************************/

	c_pt_total = cindex;

	for(index=0; index < c_pt_total; index++)
	{

	  pts[index].x = cpts[index].x;
	  pts[index].y = cpts[index].y;
	  pts[index].z = cpts[index].z;

	  point_count = c_pt_total;

	}


	printf("thru !! \n");





/************************/
/* write out new object */
/************************/

 	printf("\n\twriting out new object...");


	for(j=0; j<c_pt_total; j++)
	{
		fprintf(fout,"v %f %f %f \n",pts[j].x,
					     pts[j].y,
					     pts[j].z); 
	}



/************************/
/* close up and go home */
/************************/

	free(pts_start);
	free(texpts_start);
	free(poly_start);
	

	fclose(fin);
	fclose(fout);

	printf("done!\n");


}					/* end main */



getline(fp,buf)
FILE *fp;
char *buf;
{
	register int i,k;

again:
	i=0;

	while((k = getc(fp)) != EOF)
	{
		if(k == NEWLINE)
		{
			if(i == 0) goto again;
			*buf = 0;
			return(i);
		}
		*buf++ = k;
		++i;
	}
	return(EOF);
}




cardinal(count,tightness)
int count;
float tightness;
{

	int index;
	float t1,t2,t3;
	float f1,f2,f3,f4;
	float fx,fy,fz;
	float inc;
	int crv_index;
	int i;

sub:

	crv_index = cindex;

	for(index=1; index<=count; index++)
	{

	  t1 = index/(float)(count);
	  t2 = t1 * t1;
	  t3 = t2 * t1;

	  f1 = (0.0 - 1.0 * t1 + 2.0 * t2 - 1.0 * t3) / 2;
	  f2 = (2.0 + 0.0 * t1 - 5.0 * t2 + 3.0 * t3) / 2;
	  f3 = (0.0 + 1.0 * t1 + 4.0 * t2 - 3.0 * t3) / 2;
	  f4 = (0.0 + 0.0 * t1 - 1.0 * t2 + 1.0 * t3) / 2;

	  /* printf("f1234: %f %f %f %f\n",f1,f2,f3,f4);  */

	  fx = g[0].x * f1 + g[1].x * f2 + g[2].x * f3 + g[3].x * f4 *tightness;
	  fy = g[0].y * f1 + g[1].y * f2 + g[2].y * f3 + g[3].y * f4 *tightness;
	  fz = g[0].z * f1 + g[1].z * f2 + g[2].z * f3 + g[3].z * f4 *tightness;


	  /* printf("fxyz: %f %f %f \n",fx,fy,fz);   */


	  cpts[crv_index].x = fx;
	  cpts[crv_index].y = fy;
	  cpts[crv_index].z = fz;

	  crv_index++;
	
	}

	cindex = crv_index;

	return;

}






