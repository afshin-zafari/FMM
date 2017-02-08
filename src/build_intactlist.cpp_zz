#include "fmm.hpp"
void BuildInteractionLists_recursive(Tree &OT,int iLev,int iGroup,int jGroup){
// Recursive function
// The function evaluates the admissibility condition for a pair of groups;
// if the two groups are not 'admissible', the interactions between children
// of the two groups are checked against the admissibility condition. This
// allows not to analyze ALL combinations of group pairs (which would cost
// N^2), because groups at level iLev are considered only if the parent
// groups are neighbours).

    int nLev = length(OT) ;
    Level &ThisLev = OT.Levels(iLev) ;

    //int nGroups = length(ThisLev.group) ;
    GroupType &TestGroup   = ThisLev.group(iGroup) ;
    GroupType &SourceGroup = ThisLev.group(jGroup) ;

    double D = max(abs(SourceGroup.groupcenter-TestGroup.groupcenter)) ;    // distance (1D) between source and test groups
    if ( round(D/TestGroup.cubelength) > 1 ){
        // Low rank approximation
        OT.Levels(iLev).group(iGroup).intList.append(jGroup);
    }
    else{
        if (iLev == nLev){
            // Near field
            OT.Levels(iLev).group(iGroup).neighboursList.append(jGroup);
        }
        else{
            for ( int iSon = 1;iSon<= length(TestGroup.child);iSon++){
                for (int jSon = 1; jSon<=length(SourceGroup.child);jSon++){
                    BuildInteractionLists_recursive(OT,iLev+1,TestGroup.child(iSon),SourceGroup.child(jSon)) ;
                }
            }
        }
    }

}
void BuildInteractionLists(Tree &OT)
{
    int nLev = length(OT) ;

    elastic_vect<int> this_list ;
    this_list.append(1);
    this_list.append(1);
    for (int iLev = 1; iLev<nLev;iLev++){
        elastic_vect<int> next_list ;
        next_list.clear();
        for (int ii = 1;ii<  length(this_list);ii+=2){
            int iGroup = this_list(ii) ;
            int jGroup = this_list(ii+1) ;
            GroupType TestGroup   = OT(iLev).group(iGroup) ;
            GroupType SourceGroup = OT(iLev).group(jGroup) ;
            double D = max(abs(SourceGroup.groupcenter-TestGroup.groupcenter)) ;    // distance (1D) between source and test groups
            double r = round(D/TestGroup.cubelength);
            if ( r >= 1.0 ){
                // Low rank approximation
                OT.Levels(iLev).group(iGroup).intList.append(jGroup ) ;
            }
            else{
                if (iLev == nLev){
                    //Near field
                    OT.Levels(iLev).group(iGroup).neighboursList.append(jGroup) ;
                }
                else{
                    for (int iSon = 1;iSon < length(TestGroup.child);iSon++){
                        for (int jSon = 1;jSon< length(SourceGroup.child);jSon++){
                            next_list.append(  TestGroup.child(iSon));
                            next_list.append(SourceGroup.child(jSon));
                        }
                    }
                }
            }
        }
        this_list = next_list ;
    }

}
