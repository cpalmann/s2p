// take a series of ply files and produce a digital elevation map

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gdal.h"
#include "ogr_api.h"
#include "ogr_srs_api.h"
#include "cpl_conv.h"
#include "cpl_string.h"

#include "iio.h"
#include "lists.c"
#include "fail.c"
#include "xmalloc.c"

typedef struct {
	double x;
	double y;
	double z;
} HeightPosition;

typedef struct {
	double x;
	double y;
} Position;

struct images {
	double *cnt;
	double *pixel_value;
	Position *cell_pos;
	bool *is_empty;
	HeightPosition **heipos;
	int w, h;
};

struct ply_property {
    enum {UCHAR,FLOAT,DOUBLE,UNKNOWN} type;
    char name[0x100];
    size_t len;
};


// convert string like '28N' into a number like 32628, according to:
// WGS84 / UTM northern hemisphere: 326zz where zz is UTM zone number
// WGS84 / UTM southern hemisphere: 327zz where zz is UTM zone number
// http://www.remotesensing.org/geotiff/spec/geotiff6.html#6.3.3.1
static int get_utm_zone_index_for_geotiff(char *utm_zone)
{
    int out = 32000;
    if (utm_zone[2] == 'N')
        out += 600;
    else if (utm_zone[2] == 'S')
        out += 700;
    else
        fprintf(stderr, "error: bad utm zone value: %s\n", utm_zone);
    utm_zone[2] = '\0';
    out += atoi(utm_zone);
    return out;
}

static bool parse_property_line(struct ply_property *t, char *buf)
{
    char typename[0x100];
    bool r = 2 == sscanf(buf, "property %s %s\n", typename, t->name);
    t->type = UNKNOWN;
    if (0 == strcmp(typename, "uchar")) { t->type = UCHAR;  t->len = 1;}
    if (0 == strcmp(typename, "float")) { t->type = FLOAT;  t->len = 4;}
    if (0 == strcmp(typename, "double")){ t->type = DOUBLE; t->len = 8;}
    return r;
}

// fast forward "f" until "end_header" is found
// returns the number of 'properties'
// the array of structures *t, contains the names and sizes
// the properties in bytes, isbin is set if binary encoded
// and reads the utm zone
static size_t header_get_record_length_and_utm_zone(FILE *f_in, char *utm,
        int *isbin, struct ply_property *t)
{
    size_t n = 0;
    *isbin = 0;

    char buf[FILENAME_MAX] = {0};
    while (fgets(buf, FILENAME_MAX, f_in)) {
        if (0 == strcmp(buf, "format binary_little_endian 1.0\n")) *isbin=1;
        else if (0 == strcmp(buf, "format ascii 1.0\n")) *isbin=0;
        else {
            if (parse_property_line(t+n, buf))
                n += 1;
            else if (0 == strncmp(buf, "comment projection:", 19)) {
                sscanf(buf, "comment projection: UTM %s", utm);
            }
        }
        if (0 == strcmp(buf, "end_header\n"))
            break;
    }
    return n;
}

static void update_min_max(float *min, float *max, float x)
{
	if (x < *min) *min = x;
	if (x > *max) *max = x;
}

// re-scale a double between 0 and w-1
static int rescale_double_to_int(double x, double min, double resolution)
{
    /*double a = -(w-1)/(min-max);
    double b = min*(w-1)/(min-max);
    int r = (int) (a*x+b +0.5);*/
    int r = (int) ((x-min)/resolution+0.5);
    return r;
}

// Help to sort tabs of double
int compare (const void * a, const void * b)
{
    HeightPosition * HPa = (HeightPosition *) a;
    HeightPosition * HPb = (HeightPosition *) b;

    return HPa->z - HPb->z;
}

double weight(HeightPosition pos, Position center_pos, unsigned int flag,double pinterp)
{
    double eps=10e-3;

    double d = pow(center_pos.x-pos.x,2.0)+pow(center_pos.y-pos.y,2.0);
    d = sqrt(d);

    switch (flag)
    {
	case 6:
	    return 1.0/(d+eps);
	case 7:
	    return exp(-pinterp*d*d);
	case 8:
	    return 1/(1+exp(pinterp*d));
    }
}


// update the output images with a new height
static void add_height_to_images(struct images *x, int i, int j, Position pos, double v, int flag)
{
    uint64_t k = (uint64_t) x->w * j + i;

    switch (flag)
    {
	case -3: // just count the number of occurrences
	{
	    x->cnt[k] += 1;
	}
	break;
	case -2: // memory alloc and heights tab filling
	{
	    if (x->cnt[k])
	    {
            if (!x->heipos[k])
            {
                x->heipos[k] = xmalloc(x->cnt[k]*sizeof(HeightPosition));
                x->cnt[k]=0.0;
            }

            x->heipos[k][(int) x->cnt[k]].x = pos.x;
            x->heipos[k][(int) x->cnt[k]].y = pos.y;
            x->heipos[k][(int) x->cnt[k]].z = v;
            x->cnt[k] += 1;
	    }
	}
	break;
    }
}

static void synth_heights(struct images *x, int i, int j, int flag,
                          double pinterp, double resolution)
{
    uint64_t k = (uint64_t) x->w * j + i;

    switch (flag)
    {
        case 1: // mean
        {
            if (!x->is_empty[k])
            {
                double sum=0.,n=0.,dist;
                for(int t=0;t<x->cnt[k];t++)
                {
                    dist = MAX( abs(x->heipos[k][t].x-x->cell_pos[k].x) ,
                                abs(x->heipos[k][t].y-x->cell_pos[k].y) );
                    if (dist <= resolution/2.0)
                    {
                        sum += x->heipos[k][t].z;
                        n++;
                    }
                }
                x->pixel_value[k] = sum / n;
            }
        }
        break;
        case 2: // var
        {
            if (!x->is_empty[k])
            {
                double sum1=0.,sumC=0.;
                for(int t=0;t<x->cnt[k];t++)
                {
                    sum1 += x->heipos[k][t].z;
                    sumC += pow( x->heipos[k][t].z,2.0);
                }
                double m1 = sum1 / ( x->cnt[k] );
                double mc = sumC / ( x->cnt[k] );

                x->pixel_value[k] = mc-m1*m1;
            }
        }
        break;
        case 3: // min
        {
            if (!x->is_empty[k])
            {
                qsort (x->heipos[k], (int) x->cnt[k], sizeof(HeightPosition), compare);
                x->pixel_value[k] = x->heipos[k][0].z;
            }
        }
        break;
        case 4: // max
        {
            if (!x->is_empty[k])
            {
                qsort (x->heipos[k], (int) x->cnt[k], sizeof(HeightPosition), compare);
                x->pixel_value[k] = x->heipos[k][(int) x->cnt[k]-1].z;
            }
        }
        break;
        case 5: // median
        {
            if (!x->is_empty[k])
            {
                qsort (x->heipos[k], (int) x->cnt[k], sizeof(HeightPosition), compare);
                x->pixel_value[k] = x->heipos[k][(int) x->cnt[k]/2].z;
            }
        }
        break;
    }

    if (flag>=6)// weighted by dist from center of cell
    {
        double w;
        double sum=0.0,weighted_moy=0.0;
        int found = 0;
        if (x->is_empty[k])
        {
            for(int t=0;t<x->cnt[k];t++)
            {
              w = weight(x->heipos[k][t],x->cell_pos[k],flag,pinterp);
              sum += w;
              weighted_moy += w * x->heipos[k][t].z;
            }
            x->pixel_value[k] = weighted_moy/sum;
        }
        else
            synth_heights(x, i, j, 1, pinterp, resolution);
    }
}

int get_record(FILE *f_in, int isbin, struct ply_property *t, int n, double *data){
    int rec = 0;
    if(isbin) {
        for (int i = 0; i < n; i++) {
            switch(t[i].type) {
                case UCHAR: {
                            unsigned char X;
                            rec += fread(&X, 1, 1, f_in);
                            data[i] = X;
                            break; }
                case FLOAT: {
                            float X;
                            rec += fread(&X, sizeof(float), 1, f_in);
                            data[i] = X;
                            break; }
                case DOUBLE: {
                             double X;
                             rec += fread(&X, sizeof(double), 1, f_in);
                             data[i] = X;
                             break; }
            }
        }
    }
    else {
        int i=0;
        while (i < n && !feof(f_in)) {
            rec += fscanf(f_in,"%lf", &data[i]);  i++;
        }
    }
    return rec;
}


// open a ply file, read utm zone in the header, and update the known extrema
static void parse_ply_points_for_extrema(float *xmin, float *xmax, float *ymin,
		float *ymax, char *utm, char *fname)
{
	FILE *f = fopen(fname, "r");
	if (!f) {
		fprintf(stderr, "WARNING: can not open file \"%s\"\n", fname);
		return;
	}

	int isbin=0;
	struct ply_property t[100];
	size_t n = header_get_record_length_and_utm_zone(f, utm, &isbin, t);
	//fprintf(stderr, "%d\n", n);
	//fprintf(stderr, "%s\n", utm);

	double data[n];
	while ( n == get_record(f, isbin, t, n, data) ) {
		update_min_max(xmin, xmax, data[0]);
		update_min_max(ymin, ymax, data[1]);
	}
	fclose(f);
}

// open a ply file, and accumulate its points to the image
static void add_ply_points_to_images(struct images *img,
        double xmin, double xmax, double ymin, double ymax, int radius, double resolution,
        char utm_zone[3], char *fname, int col_idx, unsigned int flag)
{
	FILE *f = fopen(fname, "r");
	if (!f) {
		fprintf(stderr, "WARNING (from add_ply_points_to_images) : can not open file \"%s\"\n", fname);
		return;
	}

	// check that the utm zone is the same as the provided one
	char utm[3];
	int isbin=1;
	struct ply_property t[100];
	size_t n = header_get_record_length_and_utm_zone(f, utm, &isbin, t);
	if (0 != strncmp(utm_zone, utm, 3))
		fprintf(stderr, "error: different UTM zones among ply files\n");

	if (col_idx < 2 || col_idx > 5)
		exit(fprintf(stderr, "error: bad col_idx %d\n", col_idx));


	double *data = (double *) malloc(n*sizeof(double));
	double center_x,center_y,d;
	while ( n == get_record(f, isbin, t, n, data) ) {

        Position pos;
	    pos.x=data[0];
	    pos.y=-data[1];

	    double dist;
	    uint64_t k;
        int ii = rescale_double_to_int(pos.x, xmin, resolution);
	    int jj = rescale_double_to_int(pos.y, -ymax, resolution);
        int neigh = radius+5;
	    for (int j = jj-neigh; j < jj+neigh; j++)
	      for (int i = ii-neigh; i < ii+neigh; i++)
            if ( (i>=0) && (i<img->w) && (j>=0) && (j<img->h) )
        /*for (int j = 0; j < img->h; j++)
	      for (int i = 0; i < img->w; i++)*/
            {
                k = img->w * j + i;

                dist = MAX( abs(pos.x-img->cell_pos[k].x) , abs(pos.y-img->cell_pos[k].y) );

                if (dist <= resolution/2.0)
                {
                    img->is_empty[k]=false;
                    //printf("[%d %d]  A%d %dA\n",ii,jj,i,j);
                }

                if (dist <= resolution/2.0 + radius*resolution)
                {
                    //printf("(%d %d) %f %f %f %f   %f %f\n",i,j,img->cell_pos[k].x,img->cell_pos[k].y,pos.x,pos.y,dist,resolution/2.0 + radius*resolution);
                    if (col_idx == 2) {
                        add_height_to_images(img, i, j, pos, data[2], flag);
                        assert(isfinite(data[2]));
                    }
                    else
                    {
                        unsigned int rgb = data[col_idx];
                        add_height_to_images(img, i, j, pos, rgb, flag);
                    }
                }
            }
	}

    free(data);
	fclose(f);
}

void help(char *s)
{
	fprintf(stderr, "usage:\n\t"
			"ls files | %s [-c column] [-bb \"xmin xmax ymin ymax\"] [-flag flag] [-radius radius] [-param_inter param_inter] resolution out.tif\n", s);
	fprintf(stderr, "\t the resolution is in meters per pixel\n");
}

#include "pickopt.c"

int main(int c, char *v[])
{
	int col_idx = atoi(pick_option(&c, &v, "c", "2"));
	char *bbminmax = pick_option(&c, &v, "bb", "");
	int flag = atoi(pick_option(&c, &v, "flag", "1"));
	int radius = atoi(pick_option(&c, &v, "radius", "1"));
	double param_inter = atof(pick_option(&c, &v, "pinterp", "1"));

	// process input arguments
	if (c != 3) {
		help(*v);
		return 1;
	}
	float resolution = atof(v[1]);
	char *filename_out = v[2];

	// initialize x, y extrema values
	float xmin = INFINITY;
	float xmax = -INFINITY;
	float ymin = INFINITY;
	float ymax = -INFINITY;

	// process each filename from stdin to determine x, y extremas and store the
	// filenames in a list of strings, to be able to open the files again
	char fname[FILENAME_MAX], utm[3];
	struct list *l = NULL;
	while (fgets(fname, FILENAME_MAX, stdin))
	{
		strtok(fname, "\n");
		l = push(l, fname);
		parse_ply_points_for_extrema(&xmin, &xmax, &ymin, &ymax, utm, fname);
	}
	if (0 != strcmp(bbminmax, "") ) {
		sscanf(bbminmax, "%f %f %f %f", &xmin, &xmax, &ymin, &ymax);
	}
	fprintf(stderr, "xmin: %20f, xmax: %20f, ymin: %20f, ymax: %20f\n", xmin, xmax, ymin, ymax);

	double xmin_orig=xmin;
	double xmax_orig=xmax;
	double ymin_orig=ymin;
	double ymax_orig=ymax;

	// compute output image dimensions
	int w = 1 + (xmax - xmin) / resolution;
	int h = 1 + (ymax - ymin) / resolution;

	if (flag>=6) // interpolation : need more data
	{
        double mult=2.0;
        double offset = resolution/2.0 + mult*radius*resolution;
	    xmin -= offset;
	    xmax += offset;
	    ymin -= offset;
	    ymax += offset;
	    fprintf(stderr, "interpolation --> xmin: %20f, xmax: %20f, ymin: %20f, ymax: %20f\n", xmin,xmax,ymin,ymax);
	}
	else
	    radius = 0;

	struct list *begin = l;

	struct images img;
	img.w = w;
	img.h = h;
	img.cnt = xmalloc((uint64_t) w*h*sizeof(double));
	img.is_empty = xmalloc((uint64_t) w*h*sizeof(bool));
	img.pixel_value = xmalloc((uint64_t) w*h*sizeof(double));
	img.cell_pos = xmalloc((uint64_t) w*h*sizeof(Position));
	img.heipos = xmalloc((uint64_t) w*h*sizeof(HeightPosition *));

	for (int j = 0; j < img.h; j++)
	  for (int i = 0; i < img.w; i++)
	    {
		uint64_t k = img.w * j + i;
		img.cell_pos[k].x = xmin_orig + i*resolution;
		img.cell_pos[k].y = -ymax_orig + j*resolution;
		img.cnt[k] = 0;
		img.is_empty[k] = true;
		img.pixel_value[k] = 0;
		img.heipos[k] = NULL;
	    }

	printf("%d %d %d\n",w,h,img.is_empty[0]);

	l=begin;
	int n=0;
	while (l != NULL)
	{
	    n++;
	    add_ply_points_to_images(&img, xmin, xmax, ymin, ymax, radius,resolution, utm, l->current, col_idx,-3);
	    l = l->next;
	}
	l=begin;
	n=0;
	while (l != NULL)
	{
	    n++;
	    add_ply_points_to_images(&img, xmin, xmax, ymin, ymax, radius,resolution, utm, l->current, col_idx,-2);
	    l = l->next;
	}

	// heights synthesis
	for (int j = 0; j < img.h; j++)
	  for (int i = 0; i < img.w; i++)
	    synth_heights(&img,i,j,flag,param_inter,resolution);


	// set unknown values to NAN
	for (uint64_t i = 0; i < (uint64_t) img.w*img.h; i++)
		if (!img.cnt[i])
			img.pixel_value[i] = NAN;

    // open the input image
    GDALAllRegister();
    GDALDatasetH hDstDS;
    char **papszOptions = NULL;
    const char *pszFormat = "GTiff";
    GDALDriverH hDriver = GDALGetDriverByName( pszFormat );

    hDstDS = GDALCreate( hDriver, filename_out,
                         img.w, img.h, 1, GDT_Float64,
                         papszOptions );

    double adfGeoTransform[6] = { xmin_orig, resolution, 0, ymax_orig, 0, -resolution };
    OGRSpatialReferenceH hSRS;
    char *pszSRS_WKT = NULL;
    GDALRasterBandH hBand;

    GDALSetGeoTransform( hDstDS, adfGeoTransform );
    hSRS = OSRNewSpatialReference( NULL );

    char utmNumber[2];
    utmNumber[0] = utm[0];
    utmNumber[1] = utm[1];
    int nZone = atoi(utmNumber);
    int bNorth;

    if (utm[2] == 'N')
    {
        bNorth = TRUE;
    }
    else
    {
        bNorth = FALSE;
    }

    OSRSetUTM( hSRS, nZone, bNorth );
    OSRSetWellKnownGeogCS( hSRS, "WGS84" );
    OSRExportToWkt( hSRS, &pszSRS_WKT );
    OSRDestroySpatialReference( hSRS );
    GDALSetProjection( hDstDS, pszSRS_WKT );
    CPLFree( pszSRS_WKT );
    hBand = GDALGetRasterBand( hDstDS, 1 );
    GDALRasterIO( hBand, GF_Write, 0, 0, img.w, img.h,
                  img.pixel_value, img.w, img.h, GDT_Float64,
                  0, 0 );
    /* Once we're done, close properly the dataset */
    GDALClose( hDstDS );

	free(img.cnt);
	free(img.is_empty);
	free(img.pixel_value);
	free(img.cell_pos);
	for (uint64_t i = 0; i < (uint64_t) img.w*img.h; i++)
      if (img.heipos[i])
	    free(img.heipos[i]);
	free(img.heipos);

	return 0;
}
