#include "fmm.hpp"

Tree & BuildOctree(Matrix &pts, double MinCubeSide,Tree &Octree){
    /*
//BUILDOCTREE Build an Octree representation from a set of points in space.
//   OCTREE = BUILDOCTREE(X, Y, Z, MINCUBESIDE) will create a structured variable
//   containing the octree representation of (X, Y, Z).  MINCUBESIDE is the
//   smallest allowable side length of a cube in the tree.  Only non-empty
//   groups will be represented in the tree.
//
//   The structure, OCTREE, will have the following fields,
//
//   OCTREE(1:NLEVEL+1).GROUP(1:NGROUP).CHILD(1:NCHILD)
//   OCTREE(1:NLEVEL+1).GROUP(1:NGROUP).GROUPCENTER
//   OCTREE(1:NLEVEL+1).GROUP(1:NGROUP).CUBELENGTH
//   At every level, NLEVEL, of the tree, there are NGROUPS.  For each group,
//   there are NCHILD children.  GROUPCENTER is the (x0(k), y0(k), z0(k)) coordinate
//   at the center of GROUP(k). CUBELENGTH is the side length of the cube whose
//   center is at (x0(k), y0(k), z0(k)).
//
//   NOTE: A typical octree representation is numbered 0:NLEVEL. The output OCTREE in
//         BUILDOCTREE is numbered 1:NLEVEL+1 since a Matlab array can not be indexed
//         by 0.
//
//   David Sall
//   david.sall@gmail.com
//   03/01/2011
//

//////////
//////////   EQUATIONS NEEDED TO SET UP THE OCTREE
//////////

// ... The following equations* are used to generate the Octree.  At every level, except NLEVEL,
// ... of the tree the parents and children can be identified by a unique ID.  At level, NLEVEL,
// ... the children are identified by their integer value K, i.e., (x(K), y(K), z(K)) K = 1,
// ... length(x).

// ... 1) Compute the unique indicies for each element at the finest level, l = Nlevel.
// ...
// ...    i = (x - xmin) / dl
// ...    j = (y - ymin) / dl
// ...    k = (z - zmin) / dl
// ...
// ... 2) Compute a unique ID for a given (i, j, k) index at any level, ilevel.
// ...
// ...    ID = 2^(2*ilevel) * k + 2^(ilevel) * j + i
// ...
// ... 3) Compute the indicies (i, j, k) given the unique ID at any level, ilevel.
// ...
// ...    k = ID / 2^(2*ilevel)
// ...    j = ( ID - 2^(2*ilevel) * k ) / 2^(ilevel)
// ...    i = ID = 2^(2*ilevel) * k - 2^(ilevel) * j
// ...
// ... 4) Compute the PARENT (i, j, k) indicies given the CHILD (i, j, k) indicies.
// ...
// ...    i_parent = i_child / 2
// ...    j_parent = j_child / 2
// ...    k_parent = k_child / 2
// ...
// ... 5) Compute the center (x, y, z) location of a group at any level using the minimum
// ...    location of the box at the root level and the (i, j, k) indicies.
// ...
// ...    x_center = xmin + (i + .5) * dl
// ...    y_center = ymin + (j + .5) * dl
// ...    z_center = zmin + (k + .5) * dl
// ...
// ... * "Integral Equation Methods for Computational Electromagnetics", Prof. Stephen Gedney,

//////////
////////// SET UP THE BOUNDING CUBE, INITIALIZE VARIABLES, AND NUMBER OF LEVELS IN THE TREE
//////////
*/
// ... determine the point set extrema.
Matrix &xmax = max(pts,1);
Matrix &xmin = min(pts,1);

// ... expand the extrema by a small amount to make sure that, numerically, all points in the data set will
// ... be included in the tree.
Matrix &xdel = (xmax-xmin)/2*.00001;
xmax = xmax+xdel;
xmin = xmin-xdel;

// ... side length of a cube that will completely enclose all points in the data set.
double Lc = max(xmax-xmin)(1);

// ... shift the center of the bounding cube.
Matrix &XYZCenter  = (xmax+xmin)/2;
Matrix &CubeCenter = xmin+Lc/2;
xmin       = xmin-(CubeCenter-XYZCenter);

// ... number of levels in the Octree.  Level number NLEVEL is the finest level, i.e., the smallest size
// ... cube.  Level number 0 is the root level and is the largest side, i.e., bounding, cube.
// ... NOTE: Matlab does not allow idexing an array with 0 so the Octree index is shifted by 1 (so that
// ...       Octree(0:Nlevel) is indexed, Octree(1:Nlevel+1)).
double Nlevel = round(log10(Lc/MinCubeSide)/log10(2));

// ... integer constants 2^l  and 2^(2l).
double twol  = pow(2,Nlevel);
double two2l = pow(twol,2);
// ... side length of each cube at the finest level.
double dl    = Lc/twol;


//////////
////////// SET UP ARRAYS AT THE FINEST LEVEL
//////////

// ... NOTE: at the finest level, the CHILD ID's are actually the index of each point in the data set.  The
// ...       associated PARENT ID's are the ID number of the cube where each point resides.

// ... compute the parent ID for each data point ...
Matrix &i      = floor((pts.get_column(1)-xmin(1))/dl);
Matrix &j      = floor((pts.get_column(2)-xmin(2))/dl);
Matrix &k      = floor((pts.get_column(3)-xmin(3))/dl);
Matrix &parent = two2l*k+twol*j+i;
pts.print("pts");
i.print("i");
j.print("j");
k.print("k");
// ... and the center of each cube containing each point.
Matrix &groupcenter = * new Matrix (i.rows(),3);
//[xmin(1)+(i+.5)*dl xmin(2)+(j+.5)*dl xmin(3)+(k+.5)*dl];
groupcenter.set_column(1,xmin(1)+(i+.5)*dl);
groupcenter.set_column(2,xmin(2)+(j+.5)*dl);
groupcenter.set_column(3,xmin(3)+(k+.5)*dl);



//////////
////////// BUILD THE OCTREE
//////////

Matrix &child = *new Matrix (parent.rows(),parent.cols());
// Loop from the finest level to the root level.
for (int ilevel = Nlevel+1;ilevel>=1;ilevel --){

    // ... sort the PARENT ID's into ascending order.  The index sort vector are the associated CHILD
    // ... group numbers.

    ////     [parent, child] = sort(parent);
    child.rebuild(parent.rows(),parent.cols());
    sort(parent,child);
    cout << "level:" << ilevel <<endl;
    child.print();
    // ... find the index of the last parent in each group.
    ////   group = find(diff([parent; parent(end)+1]));
    #define end -1
    Matrix &group = find(diff(parent.append(parent(end)+1)));

    // ... number of groups.
    int Ngroups = length(group);

    // ... load the Octree at ILEVEL.  The PARENT ID is replaced by an integer counter (the unique ID
    // ... at each level is only used to determine the groups at each level).
    int ib = 1;
    for (int igroup = 1; igroup<=Ngroups;igroup ++){

        // ... update the ending index number of the current group.
        cout << "grp:" << igroup << "/" << Ngroups << endl;
        int ie = group(igroup);

        // ... children.
        TempIntVector rng;
        for(int i=ib;i<=ie;i++)rng.append(i);

        Octree(ilevel).group(igroup).child = child.pick_rows(rng);
        //Octree(ilevel).group(igroup).print_lists();

        // ... center of the cube containing all of the children (i.e., the cube center associated
        // ... with the PARENT ID).
        Octree(ilevel).group(igroup).groupcenter = groupcenter.get_row(child(ib));

        // ... side length of the cube whose center is given above.
        Octree(ilevel).group(igroup).cubelength  = dl;

        // ... update the starting index number of the next group.
        ib = ie+1;
    }

    // ... extract the unique PARENT ID's along with their associated indicies.
    ////[parent, igroup] = unique(parent);
    TempIntVector &igroup=*new TempIntVector;
    unique ( parent, parent,igroup);

    i = i.pick_rows(child.pick_rows(igroup));
    j = j.pick_rows(child.pick_rows(igroup));
    k = k.pick_rows(child.pick_rows(igroup));

    // ... exit the loop if at the root level.  All work is done at this point and the following
    // ... computations are unnecessary.
    if (ilevel == 1) break;

    //////////
    ////////// COMPUTE THE PARENT ID'S FOR THE NEXT LEVEL UP
    //////////

    // ... update the constants used to convert (i,j,k) indicies to the PARENT ID's.
    twol  = twol/2;
    two2l = two2l/4;
    dl = dl*2;

    // ... compute the PARENT ID's for the next level using the CHILD indicies ...
    i = floor(i/2);
    j = floor(j/2);
    k = floor(k/2);
    parent = two2l*k+twol*j+i;
    // ... and the associated center of each group.
    ////   groupcenter = [xmin(1)+(i+.5)*dl xmin(2)+(j+.5)*dl xmin(3)+(k+.5)*dl];
    groupcenter= hcat(
                      hcat( xmin(1)+(i+.5)*dl,
                            xmin(2)+(j+.5)*dl
                          ),
                      xmin(3)+(k+.5)*dl
                      );

    }
    return Octree;
}


