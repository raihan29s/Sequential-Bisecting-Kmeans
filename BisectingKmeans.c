#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>

#define DIM 8
#define DATA_NUM 120000
#define MAX_CLUSTER_NUM 16
#define KMEAN_IT_NUM 5
#define THRESH_CENTER_MOVE 1.0e-16

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


struct datapoint {
    double x[DIM];  //point coords
    int cl_id;      //cluster number point belongs to
};

struct cluster {
    double oc[DIM];
    double c[DIM];
    double sse;
    double radius;
    int point_num;
    int fp_id;
};

typedef struct datapoint dp;
typedef struct cluster cc;

/* random data generation */

void generate_random_data(dp* pdata, int dim, int data_num, cc* pcc, int max_cluster_num)
{
    int i, j, cl_id;
    srand(time(NULL));
    for (i = 0; i < DATA_NUM; i++) {
        pdata[i].cl_id = 0;
        for (j = 0; j < DIM; j++ ) {
            pdata[i].x[j] = (double)rand()/(double)RAND_MAX ;
        }
    }
    for (cl_id = 0; cl_id < max_cluster_num; cl_id++) {
        for (j = 0; j < DIM; j++ ) {
            pcc[cl_id].oc[j] = 0;
            pcc[cl_id].c[j] = 0;
            pcc[cl_id].sse = 0;
            pcc[cl_id].radius = 0;
            pcc[cl_id].fp_id = -1;
        }
    }
}

/*

calculating new centers for clusters;  returns the id of last empty/zero cluster    

*/

int find_new_cluster_centers(dp* pdata, int dim, int data_num, cc* pcc, int curr_cl_num)
{
    int i, j, cl_id;
    int cl_zero_id = -1;
    for (cl_id = 0; cl_id < curr_cl_num; cl_id++) {
        pcc[cl_id].point_num = 0;
        for (j = 0; j < DIM; j++ ) {
            pcc[cl_id].oc[j] = pcc[cl_id].c[j];
            pcc[cl_id].c[j] = 0;
            pcc[cl_id].sse = 0;
            pcc[cl_id].radius = 0;
        }
    }

    for (i = 0; i < DATA_NUM; i++) {
        cl_id = pdata[i].cl_id;
        pcc[cl_id].point_num++;
        for (j = 0; j < DIM; j++ ) {
            pcc[cl_id].c[j] += pdata[i].x[j];
        }
    }
    for (cl_id = 0; cl_id < curr_cl_num; cl_id++)
      {
        if (pcc[cl_id].point_num == 0) 
	  {
            cl_zero_id = cl_id;
            continue;
	  }
        for (j = 0; j < DIM; j++ ) {
            pcc[cl_id].c[j] /= pcc[cl_id].point_num;
        }
    }
    return cl_zero_id;
    
}

/* distance between two points */

double sqr_euc_dist(double* p1, double* p2, int dim)
{
    int i;
    double dist = 0.0;
    for(i = 0; i < dim; i++) {
        dist += (p1[i] - p2[i]) * (p1[i] - p2[i]);
    }
    return dist;
}

/*
  computes SSE and radious^2 of the cluster
*/

void  compute_sse_radius_for_clusters(dp* pdata, int dim, int data_num, cc* pcc, int curr_cl_num)
{
    double dist_p_c;
    int i, cl_id;

    for (cl_id = 0; cl_id < curr_cl_num; cl_id++) {
        pcc[cl_id].sse = 0;
        pcc[cl_id].radius = 0;
        pcc[cl_id].fp_id = -1;
    }
    
    for (i = 0; i < data_num; i++) {
        cl_id = pdata[i].cl_id;
        dist_p_c = sqr_euc_dist( pdata[i].x, pcc[cl_id].c, dim);
        pcc[cl_id].sse += dist_p_c;
        if (dist_p_c < pcc[cl_id].radius)
            continue;
        pcc[cl_id].radius = dist_p_c;
        pcc[cl_id].fp_id = i;
    }
}

/*
it excludes point which is used for creation of new cluster. Also it shift center of old cluster to position 2C-p  
*/

int compute_sse_radius_for_symmetric_cluster_with_excl(dp* pdata, int dim, int data_num, cc* pcc, int cl_id, int excl_id)
{
    double dist_p_c;
    int i, j;    

    pcc[cl_id].sse = 0.0;
    pcc[cl_id].radius = 0.0;
    pcc[cl_id].point_num--;
    for (j = 0; j < dim; j++) {
        pcc[cl_id].c[j] = 2 * pcc[cl_id].c[j] - pdata[excl_id].x[j];
    }

    for (i = 0; i < data_num; i++) {
        if (pdata[i].cl_id != cl_id)
            continue;
        dist_p_c = sqr_euc_dist( pdata[i].x, pcc[cl_id].c, dim);
        pcc[cl_id].sse += dist_p_c;
        if (dist_p_c < pcc[cl_id].radius)
            continue;
        pcc[cl_id].radius = dist_p_c;
        pcc[cl_id].fp_id = i;
    }
}

/* bisect the worst cluster with worst SSE */

void bisect_worst_cluster(dp* pdata, int dim, int data_num, cc* pcc, int curr_cl_num, int new_cl_id)
{
    int worst_cl_id;
    int worst_sse = 0.0;
    int cl_id;
    int j;

    compute_sse_radius_for_clusters(pdata, dim, data_num, pcc, curr_cl_num);
    for(cl_id = 0; cl_id < curr_cl_num; cl_id++) {
        if(worst_sse > pcc[cl_id].sse)
            continue;
        worst_sse = pcc[cl_id].sse;
        worst_cl_id = cl_id;
    }

    //bisection itself
    
    //create new cluster or redefine zero one
    pcc[new_cl_id].point_num = 1;
    pcc[new_cl_id].sse = 0.0; 
    pcc[new_cl_id].radius = 0.0;
    pcc[new_cl_id].fp_id = pcc[worst_cl_id].fp_id;
    for(j = 0; j < dim; j++) {
        pcc[new_cl_id].oc[j] = 0;
        pcc[new_cl_id].c[j] = pdata[pcc[new_cl_id].fp_id].x[j];
    }
    pdata[pcc[new_cl_id].fp_id].cl_id = new_cl_id;
    //fix worst cluster
    compute_sse_radius_for_symmetric_cluster_with_excl(pdata, dim, data_num, pcc, worst_cl_id, pcc[new_cl_id].fp_id);

}

/*
  rearrange_points between cluster centers. Checking for min distance
*/

void rearrange_points(dp* pdata, int dim, int data_num, cc* pcc, int curr_cl_num)
{
    int i, cl_id, best_cl_id;
    double dist_p_c;
    double min_dist;
    #
    for (i = 0; i < data_num; i++) {
        best_cl_id = 0;
        min_dist = sqr_euc_dist(pcc[best_cl_id].c, pdata[i].x, dim);
        for(cl_id = 1; cl_id < curr_cl_num; cl_id++ ) {
            dist_p_c = sqr_euc_dist(pcc[cl_id].c, pdata[i].x, dim);
            if (dist_p_c < min_dist) {
                min_dist = dist_p_c;
                best_cl_id = cl_id;
            }
        }
        pdata[i].cl_id = best_cl_id;
    }
}

int main()
{
    dp *p;
    cc *pcc;
    int curr_cl_num, cl_id, j;
    int fix_cluster_id;
    int ret;
    int it, stop_it;
    double center_dev;

    p = (dp*) malloc(DATA_NUM * (sizeof(dp)));
    pcc = (cc*) malloc(MAX_CLUSTER_NUM * (sizeof(cc)));

    printf("Generating %d random data points started.\n", DATA_NUM);
    generate_random_data(p, DIM, DATA_NUM, pcc, MAX_CLUSTER_NUM);
    printf("Generating random data points ended.\n");

    curr_cl_num = 1;
    while (curr_cl_num++ < MAX_CLUSTER_NUM) {
        bisect_worst_cluster(p, DIM, DATA_NUM, pcc, curr_cl_num - 1 , curr_cl_num);
        it = 0;
        stop_it = 0;
        while(it < KMEAN_IT_NUM && (! stop_it)) {
            rearrange_points(p, DIM, DATA_NUM, pcc, curr_cl_num);
            ret = find_new_cluster_centers(p, DIM, DATA_NUM, pcc, curr_cl_num);
            fix_cluster_id = ret;
            if ( ret != -1) { // empty cluster found
                bisect_worst_cluster(p, DIM, DATA_NUM, pcc, curr_cl_num, fix_cluster_id);
                it = 0;
                continue;       
            }
            it++;
            center_dev = 0;
            for(cl_id = 0; cl_id < curr_cl_num; cl_id++) {
                center_dev += sqr_euc_dist(pcc[cl_id].c, pcc[cl_id].oc, DIM);
            }
            if (center_dev < THRESH_CENTER_MOVE) {
                stop_it = 1;
            }
        }
    }
     
    compute_sse_radius_for_clusters(p, DIM, DATA_NUM, pcc, curr_cl_num);

    //output: cluster_centers and cluster_size
    for (cl_id = 0; cl_id < MAX_CLUSTER_NUM; cl_id++) {
        printf("====\n");
        printf("Cluster_id = %d , fp_id = %d, cluster_radius %lf\n", cl_id, pcc[cl_id].fp_id, sqrt(pcc[cl_id].radius));
        for(j = 0; j < DIM; j++) {
            printf(" x[%d] = %lf", j, pcc[cl_id].c[j]);
            if ((j + 1) % 4) 
                continue;
            printf("\n");
        }
    }
    

	return 0;
}


