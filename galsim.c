/*
 * File: galsim.c
 * --------------
 * Barnes hut algorithm, simulation of galaxy
 *
 * HPC Assignment 1
 * Carl Christian Kirchman
 * Alessandro Piccolo
 * Anton Sjöberg
 */

#include <stdio.h>
#include <stdlib.h>
#include "graphics.h"
#include "math.h"
#include "file_operations.h"
#include <sys/time.h>

typedef struct {
  double mass;
  double x_coord;
  double y_coord;
  double x_velocity;
  double y_velocity;
} particle;

typedef struct {
  double mass;
  double x_coord;
  double y_coord;
} returnV;

typedef struct {
  double fx;
  double fy;
} force_struct;

struct treeNode_ {
  struct treeNode_* child_1;        /* Upper left                             */
  struct treeNode_* child_2;        /* Upper right                            */
  struct treeNode_* child_3;        /* Down left                              */
  struct treeNode_* child_4;        /* Down right                             */
  particle* part_i;					        /* Particle in treenode 				          */
  unsigned int leaf;				        /* 1 if leaf, 0 otherwise				          */
  double llc[2];                    /* Coordinates of lower left corner (llc) */
  int lvl;                          /* Tree depth (root is lvl 1)             */
  double CoMX;                      /* Centre of mass X coordinate            */
  double CoMY;                      /* Centre of mass Y coordinate            */
  double nodeMass;                  /* Total mass of that square              */
  /* Number of particles in subtree (stops counting at 2) */ 
  int particle_sum;	            
};
typedef struct treeNode_ treeNode;

/* Creates 4 children with: children set to NULL, part_i set to NULL */
void QuadtreeBuild(treeNode* node, particle* particles, int N);
void search_cleanup(treeNode* root);
void insert_tree(treeNode* node, particle* particle);
void create_branch(treeNode* node); 
returnV centre_of_mass(treeNode* node);
force_struct tree_force_calculation(particle* particle, treeNode* node,
                                    double G, double theta_max, double epsilon);
void print_particles(particle* particles, int N);
void print_quadTree(treeNode* rootNode);
void free_tree(treeNode* rootNode);
static double get_wall_seconds();

/* Print */
/*
const double circleRadius = 0.0015, circleColor = 0;
const int windowWidth = 1200;
*/

int main(int argc, char *argv[]) {
  if (argc != 6) { /* Five input arguments check */
    printf("Should be,\n./galsim N input_data/filename");
    printf(" nsteps delta_t theta_max\n");
    return 1;
  }
  int i, j, k;                              /* Looping variables              */
  double *p;                                /* Buffer with particle info      */
  int    const N          = atoi(argv[1]);  /* Number of particles            */
  int    const nsteps     = atoi(argv[3]);  /* Number of steps                */
  double const dt        = atof(argv[4]);   /* Time step delta t              */
  double const theta_max = atof(argv[5]);   /* Theta max threshold            */
  double const G         = 100/(double)N;   /* Gravitational constant         */
  double a[2];                              /* Acceleration of particle       */        
  double epsilon = 0.001;
  double const L = 1, W = 1;
  treeNode* rootNode;                       /* Starting root tree             */  
 
  p = (double *)malloc(N*5*sizeof(double)); /* Malloc to save stacken   	    */
  read_doubles_from_file(N*5, p, argv[2]);  /* Read particle info into buffer */
  
  /* Print */ 
  /*
  InitializeGraphics(argv[0], windowWidth, windowWidth);
  SetCAxes(0, 1); 
  */
  
  /* Get particle data into an array of particles struct */
  particle* particles = (particle*)malloc(N*sizeof(particle));
  int index_p;

  double startTime = get_wall_seconds();
  for (i = 0; i < N; i++) {                 /* Can probably unroll this loop  */
    index_p                 = i*5;
    particles[i].mass       = p[index_p];
    particles[i].x_coord    = p[index_p+1];
    particles[i].y_coord    = p[index_p+2];
    particles[i].x_velocity = p[index_p+3];
    particles[i].y_velocity = p[index_p+4];
  }

  double treebuildtime      = 0;
  double centreOfMasstime   = 0;
  double forceTime          = 0;
  double velTime            = 0;
  double coordTime          = 0;
  double startBuild, massStart, forceStart, velStart, coordStart;

  /* Time steps */
  for (k = 0; k < nsteps; k++) { 
    rootNode = (treeNode*)malloc(sizeof(treeNode));
    /* Fill tree with respective particles */
    startBuild = get_wall_seconds();
    QuadtreeBuild(rootNode, particles, N);
    treebuildtime += (get_wall_seconds()- startBuild);
    /* Calculate centre of mass and total mass for all the nodes */
    massStart = get_wall_seconds();
    centre_of_mass(rootNode);
    centreOfMasstime += (get_wall_seconds()- massStart);
    
    /* Calculate Force */
    force_struct F;
    double  mi;                             /* Mass of particle i             */
    for (i = 0; i < N; i++) {
      forceStart = get_wall_seconds();
      F = tree_force_calculation(&particles[i],rootNode,G,theta_max,epsilon);
      forceTime += (get_wall_seconds()- forceStart);
      coordStart = get_wall_seconds();
      mi   = particles[i].mass;
      a[0] = F.fx/mi;                      /* Acceleration x                  */
      a[1] = F.fy/mi;                      /* Acceleration y                  */
      particles[i].x_velocity += dt*a[0];  /* Update velocity x               */
      particles[i].y_velocity += dt*a[1];  /* Update velocity y               */
      coordTime = (get_wall_seconds() - coordStart);
    }
    
    /* Update position for each particle */
    velStart = get_wall_seconds();
    for (i = 0; i < N; i++) {
      particles[i].x_coord += dt*particles[i].x_velocity; 
      particles[i].y_coord += dt*particles[i].y_velocity;     
    }
    velTime = (get_wall_seconds() - velStart);
    
    /* Print */
    /*
    if (k%10 == 0) {
     ClearScreen();  
     for(i = 0; i < N; i++) { // Draw circles
        DrawCircle(particles[i].x_coord ,particles[i].y_coord,
                   L, W, circleRadius, circleColor);
      }
      Refresh();  
      usleep(1000); // Sleep a short while to avoid screen flickering
    } 
    */
    free_tree(rootNode);
  }

  /* Make the returnvector*/
  for(i = 0; i < N; i++) {
    p[i*5] = particles[i].mass;
    p[i*5+1] = particles[i].x_coord;
    p[i*5+2] = particles[i].y_coord;
    p[i*5+3] = particles[i].x_velocity;
    p[i*5+4] = particles[i].y_velocity;
  }

  double endTime = get_wall_seconds() - startTime;
  printf("------------------------------------------------------------\n\n");
  printf("treebuildtime    = %f \n", treebuildtime);
  printf("centreOfMasstime = %f \n", centreOfMasstime);
  printf("forceTime        = %f \n", forceTime);
  printf("velTime          = %f \n", velTime);
  printf("coordTime        = %f \n", coordTime);
  printf("\nTotal time = %f \n\n", endTime);
  printf("------------------------------------------------------------\n");

  /* Output particle info into gal file */
  write_doubles_to_file(N*5, p, "BHOutput.gal");

  /* Print */
  /*
  FlushDisplay();
  CloseDisplay();
  */

  /* Free malloc */
  free(p); 
  free(particles);
  return 0;
}

/* Free quadtree */
void free_tree(treeNode* rootNode) {
  if (rootNode->child_1) {
    free_tree(rootNode->child_1);
  }
  if (rootNode->child_2) {
    free_tree(rootNode->child_2);
  }
  if (rootNode->child_3) {
    free_tree(rootNode->child_3);
  }
  if (rootNode->child_4) {
    free_tree(rootNode->child_4);
  }
  rootNode->part_i = NULL;
  rootNode->leaf = 0;
  rootNode->lvl = 0;
  rootNode->nodeMass = 0;
  rootNode->child_1 = NULL;
  rootNode->child_2 = NULL;
  rootNode->child_3 = NULL;
  rootNode->child_4 = NULL;
  free(rootNode);
  rootNode = NULL;
}

/* Print quadtree */
void print_quadTree(treeNode* node) {
  printf("Level: %d\n", node->lvl);
  if (node->part_i != NULL) {
    printf("x = %f, y = %f\n", node->part_i->x_coord,node->part_i->y_coord);
  }
  if (node->child_1 != NULL) {
    print_quadTree(node->child_1);
  }
  if (node->child_2 != NULL) {
    print_quadTree(node->child_2);
  }
  if (node->child_3 != NULL) {
    print_quadTree(node->child_3);
  }
  if (node->child_4 != NULL) {
    print_quadTree(node->child_4);
  }
}

/* Print particles coordinates */
void print_particles(particle* particles, int N) {
  int i;
  for (i = 0; i < N; i++) {
    printf("Particle %d x_coord = %f, y_coord = %f, mass = %f\n", 
        i, particles[i].x_coord, particles[i].y_coord,particles[i].mass);
  }
}

/* Build the tree */
void QuadtreeBuild(treeNode* rootNode, particle* particles, int N) {
  /* Filling in root tree data */
  rootNode->child_1  = NULL;
  rootNode->child_2  = NULL;
  rootNode->child_3  = NULL;
  rootNode->child_4  = NULL;
  rootNode->part_i   = NULL;

  rootNode->llc[0]  = 0;
  rootNode->llc[1]  = 0;
  rootNode->lvl     = 1;
  rootNode->leaf	  = 1;
  rootNode->particle_sum = 0;
  int i;
  for (i = 0; i < N; i++){
    insert_tree(rootNode, &particles[i]);
  }
  /* Traverse and clean up nodes with 0 particles */ 
  search_cleanup(rootNode);
}

/* Insert the tree */
void insert_tree(treeNode* node, particle* particle) {
  double xCoord = particle->x_coord; 
  double yCoord = particle->y_coord;
  double shift_lvl  = (1 << node->lvl);
  double offset     = 1/shift_lvl;
  double mid_x = node->llc[0] + offset;       /* x coordinate of the mid line */
  double mid_y = node->llc[1] + offset;       /* y coordinate of the mid line */

  /* Node subtree contains more than 1 particle */
  if (node->particle_sum == 2) {   
    if (xCoord > mid_x) {                     /* Right side                   */
      if (yCoord >= mid_y) {                  /* Upper right                  */
  		insert_tree(node->child_2,particle);
      } 
      else {                                  /* Down right                   */  
        insert_tree(node->child_4,particle);
      }
    } 
    else {                                    /* Left side                    */
      if (yCoord >= mid_y) {                  /* Upper left                   */
        insert_tree(node->child_1,particle); 
      } 
      else {                                  /* Down left                    */  
        insert_tree(node->child_3,particle);
      }
    }
  }
  else if (node->particle_sum == 1) {
    node->particle_sum = 2;
    create_branch(node);
    node->leaf = 0;
    double xCoordP = node->part_i->x_coord;
    double yCoordP = node->part_i->y_coord;

    /* Moving the particle already in the node */
    if (xCoordP > mid_x) {                     /* Right side                  */
      if (yCoordP >= mid_y) {                  /* Upper right                 */
        node->child_2->part_i = node->part_i;
        node->part_i = NULL;
	node->child_2->particle_sum = 1;
      } 
      else {                                   /* Down right                  */  
        node->child_4->part_i = node->part_i;
        node->part_i = NULL;
	node->child_4->particle_sum = 1;
      }
    } 
    else {                                      /* Left side                  */
      if (yCoordP >= mid_y) {                   /* Upper left                 */
        node->child_1->part_i = node->part_i;
        node->part_i = NULL;
	node->child_1->particle_sum = 1;
    } 
      else {                                    /* Down left                  */  
        node->child_3->part_i = node->part_i;
        node->part_i = NULL;
	node->child_3->particle_sum = 1;
      }
    }
   
    /* Inserting the new particle*/
    if(xCoord > mid_x) {                         /* Right side                */
      if(yCoord >= mid_y) {                      /* Upper right               */
        insert_tree(node->child_2,particle);
      } 
      else {                                     /* Down right                */  
        insert_tree(node->child_4,particle);
      }
    } 
    else {                                       /* Left side                 */
      if(yCoord >= mid_y) {                      /* Upper left                */
        insert_tree(node->child_1,particle); 
      } 
      else {                                     /* Down left                 */  
        insert_tree(node->child_3,particle);
      }
    }
  }
  /* Subtree rooted at node is empty */ 
  else if (node->particle_sum == 0) { 
    node->part_i = particle;
    node->leaf   = 1;
    node->particle_sum = 1;
  } 
  else {
    printf("Invalid output from tree_particle_cntr\n");
  } 
  return; 
}

/* Create a branch */
void create_branch(treeNode* node) {
  double shift_lvl  = (1 << node->lvl);
  double offset     = 1/shift_lvl;
  int lvl           = node->lvl + 1;

  /* Child 1 */
  treeNode* node_ul = (treeNode*)malloc(sizeof(treeNode));
  node->child_1     = node_ul;
  node_ul->lvl      = lvl;
  node_ul->child_1  = NULL;
  node_ul->child_2  = NULL;
  node_ul->child_3  = NULL;
  node_ul->child_4  = NULL;
  node_ul->part_i   = NULL;
  node_ul->leaf     = 1;
  node_ul->particle_sum = 0;

  /* Child 2 */
  treeNode* node_ur = (treeNode*)malloc(sizeof(treeNode));
  node->child_2     = node_ur;
  node_ur->lvl      = lvl;
  node_ur->child_1  = NULL;
  node_ur->child_2  = NULL;
  node_ur->child_3  = NULL;
  node_ur->child_4  = NULL;
  node_ur->part_i   = NULL;
  node_ur->leaf     = 1;
  node_ur->particle_sum = 0;

  /* Child 3 */
  treeNode* node_dl = (treeNode*)malloc(sizeof(treeNode));
  node->child_3     = node_dl;
  node_dl->lvl      = lvl;
  node_dl->child_1  = NULL;
  node_dl->child_2  = NULL;
  node_dl->child_3  = NULL;
  node_dl->child_4  = NULL;
  node_dl->part_i   = NULL;
  node_dl->leaf     = 1;
  node_dl->particle_sum= 0;

  /* Child 4 */
  treeNode* node_dr = (treeNode*)malloc(sizeof(treeNode));
  node->child_4     = node_dr;
  node_dr->lvl      = lvl;
  node_dr->child_1  = NULL;
  node_dr->child_2  = NULL;
  node_dr->child_3  = NULL;
  node_dr->child_4  = NULL;
  node_dr->part_i   = NULL;
  node_dr->leaf     = 1;
  node_dr->particle_sum = 0;

  /* Calculating the lower left cornet coordinates for the children */
  double llc_0      = node->llc[0];
  double llc_1      = node->llc[1];
  node_ul->llc[0]   = llc_0;
  node_ul->llc[1]   = llc_1 + offset; 
  node_ur->llc[0]   = llc_0 + offset; 
  node_ur->llc[1]   = llc_1 + offset;
  node_dl->llc[0]   = llc_0;
  node_dl->llc[1]   = llc_1;  
  node_dr->llc[0]   = llc_0 + offset;
  node_dr->llc[1]   = llc_1;
}

/* Clean up the empty leaves */
void search_cleanup(treeNode* node) {
  /* Traverse the tree */
  if (node->child_1->leaf != 1) {
    search_cleanup(node->child_1);
  }
  if (node->child_2->leaf != 1) {
    search_cleanup(node->child_2);
  }
  if (node->child_3->leaf != 1) {
    search_cleanup(node->child_3);
  }
  if (node->child_4->leaf != 1) {
    search_cleanup(node->child_4);
  }

  /* Remove the empty leaves */ 
  if (node->child_1->leaf == 1 && node->child_1->part_i==NULL) {
    free(node->child_1);  /* Free memory from heap 	    */
    node->child_1 = NULL;   /* Represents a removed child */
  }
  if(node->child_2->leaf == 1 && node->child_2->part_i==NULL) {
    free(node->child_2);
    node->child_2 = NULL;
  }
  if(node->child_3->leaf == 1 && node->child_3->part_i==NULL) {
    free(node->child_3);
    node->child_3 = NULL;
  }
  if(node->child_4->leaf == 1 && node->child_4->part_i==NULL) {
    free(node->child_4);
    node->child_4 = NULL;
  } 
}

/* Calculate centre of mass and total mass for each node */
returnV centre_of_mass(treeNode* node) { /* Returns: [mass,Xcoord,Ycoord] */
  returnV returnVector;
  if (node->particle_sum == 1) {
    node->nodeMass = node->part_i->mass;
    node->CoMX = node->part_i->x_coord;
    node->CoMY = node->part_i->y_coord;
    returnVector.mass = node->nodeMass;
    returnVector.x_coord = node->CoMX;
    returnVector.y_coord = node->CoMY;
    return returnVector;
  }
  else {
    double mass_sum = 0.0, CoMY_sum = 0.0, CoMX_sum = 0.0;
    if (node->child_1 != NULL) {
      returnV data_child_1 = centre_of_mass(node->child_1);
      mass_sum += data_child_1.mass;
      CoMX_sum += data_child_1.x_coord*data_child_1.mass;
      CoMY_sum += data_child_1.y_coord*data_child_1.mass;
    }
    if (node->child_2 != NULL) {
      returnV data_child_2 = centre_of_mass(node->child_2);
      mass_sum += data_child_2.mass;
      CoMX_sum += data_child_2.x_coord*data_child_2.mass;
      CoMY_sum += data_child_2.y_coord*data_child_2.mass;      
    }
    if (node->child_3 != NULL) {
      returnV data_child_3 = centre_of_mass(node->child_3);
      mass_sum += data_child_3.mass;
      CoMX_sum += data_child_3.x_coord*data_child_3.mass;
      CoMY_sum += data_child_3.y_coord*data_child_3.mass;
    }
    if (node->child_4 != NULL) {
      returnV data_child_4 = centre_of_mass(node->child_4);
      mass_sum += data_child_4.mass;
      CoMX_sum += data_child_4.x_coord*data_child_4.mass;
      CoMY_sum += data_child_4.y_coord*data_child_4.mass;
    } 
    /* Store info in node */    
    node->nodeMass = mass_sum;
    node->CoMX = CoMX_sum/mass_sum;
    node->CoMY = CoMY_sum/mass_sum;
    /* Update structure for return value */   
    returnVector.mass    = node->nodeMass;
    returnVector.x_coord = node->CoMX;
    returnVector.y_coord = node->CoMY;

    return returnVector;
  }
}

/* Calculate force for particle from particles in node */
force_struct tree_force_calculation(particle* particle, treeNode* node,
                                    double G, double theta_max, double epsilon){ 
  force_struct force_container;
  force_container.fx = 0; force_container.fy = 0; /* Needs to be init */
  double m = particle->mass;
  double r_bold[2], r;
  if (node->part_i != NULL) { 
  /* Can only be one particle in subtree if the node has one particle*/
	  if (particle != node->part_i) { /* Should not interact with itself */
      r_bold[0] =  particle->x_coord - node->CoMX;
	    r_bold[1] =  particle->y_coord - node->CoMY;
      r         = sqrt( r_bold[0]*r_bold[0] + r_bold[1]*r_bold[1] );
      /* Always calculate F with Plummer spheres */
      force_container.fx -=  G*m*node->nodeMass*r_bold[0]/
      ((r+epsilon)*(r+epsilon)*(r+epsilon));
      force_container.fy -=  G*m*node->nodeMass*r_bold[1]/
      ((r+epsilon)*(r+epsilon)*(r+epsilon));
    }
  }
  else {
    r_bold[0]   = particle->x_coord - node->CoMX;
	  r_bold[1]   = particle->y_coord - node->CoMY;
    r           = sqrt( r_bold[0]*r_bold[0] + r_bold[1]*r_bold[1] );
    int sublvl  = node->lvl - 1;
    double shifted_sublvl = 1<<sublvl;
    double D    = 1/shifted_sublvl; /* Size of box node is in*/
    if (D/r < theta_max){
      /* Always calculate F with Plummer spheres */
      force_container.fx -=  G*m*node->nodeMass*r_bold[0]/
      ((r+epsilon)*(r+epsilon)*(r+epsilon));
      force_container.fy -=  G*m*node->nodeMass*r_bold[1]/
      ((r+epsilon)*(r+epsilon)*(r+epsilon));
    } 
    else {
      if (node->child_1!=NULL) {
        force_struct struct1 = tree_force_calculation(particle,node->child_1,G,
                                                      theta_max,epsilon);
        force_container.fx = force_container.fx + struct1.fx;
        force_container.fy = force_container.fy + struct1.fy; 
      }
      if (node->child_2!=NULL) {
        force_struct struct2 = tree_force_calculation(particle,node->child_2,G,
                                                      theta_max,epsilon);
        force_container.fx = force_container.fx + struct2.fx;
        force_container.fy = force_container.fy + struct2.fy; 
      }
      if (node->child_3!=NULL) {
        force_struct struct3 = tree_force_calculation(particle,node->child_3,G,
                                                      theta_max,epsilon);
        force_container.fx = force_container.fx + struct3.fx;
        force_container.fy = force_container.fy + struct3.fy; 
      }
      if (node->child_4!=NULL) {
        force_struct struct4 = tree_force_calculation(particle,node->child_4,G,
                                                      theta_max,epsilon);
        force_container.fx = force_container.fx + struct4.fx;
        force_container.fy = force_container.fy + struct4.fy; 
      }
    }
  }
  return force_container;
}

static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}
