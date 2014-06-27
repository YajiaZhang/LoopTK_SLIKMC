/*
    LoopTK: Protein Loop Kinematic Toolkit
    Copyright (C) 2007 Stanford University

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

//!----------------------------------------------------------------------
//! Some of this code was adapted from code which is protected by the following copyright:
//! 	Copyright (C) 2003 
//!      Chaok Seok, Evangelos Coutsias, Matthew Jacobson, and Ken Dill.
//!      UCSF and Univeristy of New Mexico.
//! 	Witten by Chaok Seok 2003.  
//! This does not prevent distribution of LoopTK.
//!----------------------------------------------------------------------
//!-----------------------------------------------------------------------
#include "PIKAlgorithms.h"
#include "PConstants.h"
#include "PExtension.h"
#include "PSturm.h"
#include "PTripepClosure.h"
#include "PTools.h"
#define PI 3.14159265

IKSolutions PExactIKSolver::FindSolutions(PProtein *loop,int DOF_indices_to_use[3]) {
  Vector3 endPriorG;
  Vector3 endG;
  Vector3 endNextG;
  endPriorG = loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[2])->getPos();
  endG = loop->getAtomAtRes(PID::C,DOF_indices_to_use[2])->getPos();
  endNextG = loop->getAtomAtRes(PID::O,DOF_indices_to_use[2])->getPos();
  return FindSolutions(loop,DOF_indices_to_use,&endPriorG,&endG,&endNextG);
}

//Yajia added this function.
int PExactIKSolver::FindSolutions(PProtein *loop, int DOF_indices_to_use[3], Vector3 *endPriorG, Vector3 *endG, Vector3 *endNextG, IKSolutions& solutions){
//!-----------------------------------------------------------------------
//! This is a sample driver routine to reconstruct tripeptide loops
//! from the coordinates in a pdb file using a canonical bond lengths and
//! angles.
//!-----------------------------------------------------------------------
//  use in_out
//  use tripep_closure
//!-----------------------------------------------------------------------
//  implicit none
//  real(dp) :: time1, time2
//  character(len=100) :: in_pdb, out_prefix, out_pdb
  char in_pdb[100], out_prefix[100], out_pdb[100];
//  character(len=4) :: res_name(5)
  char res_name[5][4];
//  integer :: n0, n_soln, unit = 17, i, j, k, n1, n2, order(max_soln)
  int n0, n_soln, i, j, k, n1, n2;
//  real(dp) :: r_n(3,5), r_a(3,5), r_c(3,5)
  double r_n[5][3], r_a[5][3], r_c[5][3];
//  real(dp) :: r_soln_n(3,3,max_soln), r_soln_a(3,3,max_soln), r_soln_c(3,3,max_soln)
  double r_soln_n[max_soln][3][3], r_soln_a[max_soln][3][3], r_soln_c[max_soln][3][3];
//  real(dp) :: rmsd, sum, dr(3)
  double rmsd, sum, dr[3];
//  real(dp) :: r0_n(3,3), r0_a(3,3), r0_c(3,3)
  double r0_n[3][3], r0_a[3][3], r0_c[3][3];
  
//  logical :: calc_rmsd, write_out_pdb
  int  calc_rmsd, write_out_pdb;
//  real(dp) :: b_len(6), b_ang(7), t_ang(2)
  double b_len[6], b_ang[7], t_ang[2];

//  ! input parameters (lengths)
//  real(dp), parameter :: b_na = 1.45d0, b_ac = 1.52d0, b_cn = 1.33d0
  double b_na = 1.45e0, b_ac = 1.52e0, b_cn = 1.33e0;
//  ! input parameters (angles in radians)
//  real(dp), parameter :: ang_nac = 111.6*deg2rad, ang_acn = 117.5*deg2rad
  double ang_nac = 112*deg2rad, ang_acn = 120.0*deg2rad;
//  real(dp), parameter :: ang_cna = 120.0*deg2rad
  double ang_cna = 121.9*deg2rad;

//!-----------------------------------------------------------------------

//  call my_timer(time1)

//  !! initialize: set up input parameters
//  ! bond lengths
//  b_len(1:6) = (/ b_ac, b_cn, b_na, b_ac, b_cn, b_na /)
  //cout<<"InFindSolutions2"<<endl;
  Vector3 p1,p2;
  p1=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[0])->getPos();
  p2=loop->getAtomAtRes(PID::C,DOF_indices_to_use[0])->getPos();
  b_len[0]=p1.distance(p2);
  
  p1=loop->getAtomAtRes(PID::C,DOF_indices_to_use[0])->getPos();
  p2=loop->getAtomAtRes(PID::N,DOF_indices_to_use[1])->getPos();
  b_len[1]=p1.distance(p2);
  p1=loop->getAtomAtRes(PID::N,DOF_indices_to_use[1])->getPos();
  p2=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[1])->getPos();
  b_len[2]=p1.distance(p2);
  
  p1=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[1])->getPos();
  p2=loop->getAtomAtRes(PID::C,DOF_indices_to_use[1])->getPos();
  b_len[3]=p1.distance(p2);
  
  p1=loop->getAtomAtRes(PID::C,DOF_indices_to_use[1])->getPos();
  p2=loop->getAtomAtRes(PID::N,DOF_indices_to_use[2])->getPos();
  b_len[4]=p1.distance(p2);
  
  p1=loop->getAtomAtRes(PID::N,DOF_indices_to_use[2])->getPos();
  p2=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[2])->getPos();
  b_len[5]=p1.distance(p2);
  
  
//  ! bond angles
//  b_ang(1:7) = (/ ang_nac, ang_acn, ang_cna, ang_nac, ang_acn, ang_cna, ang_nac /)

  Vector3 q;
  p1=loop->getAtomAtRes(PID::N,DOF_indices_to_use[0])->getPos();
  p2=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[0])->getPos();
  q=loop->getAtomAtRes(PID::C,DOF_indices_to_use[0])->getPos();
  b_ang[0]=AngleBetweenVectors(p2-p1,p2-q);
  
  p1=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[0])->getPos();
  p2=loop->getAtomAtRes(PID::C,DOF_indices_to_use[0])->getPos();
  q=loop->getAtomAtRes(PID::N,DOF_indices_to_use[1])->getPos();
  b_ang[1]=AngleBetweenVectors(p2-p1,p2-q);
  
  p1=loop->getAtomAtRes(PID::C,DOF_indices_to_use[0])->getPos();
  p2=loop->getAtomAtRes(PID::N,DOF_indices_to_use[1])->getPos();
  q=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[1])->getPos();
  b_ang[2]=AngleBetweenVectors(p2-p1,p2-q);
  
  p1=loop->getAtomAtRes(PID::N,DOF_indices_to_use[1])->getPos();
  p2=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[1])->getPos();
  q=loop->getAtomAtRes(PID::C,DOF_indices_to_use[1])->getPos();
  b_ang[3]=AngleBetweenVectors(p2-p1,p2-q);
  
  p1=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[1])->getPos();
  p2=loop->getAtomAtRes(PID::C,DOF_indices_to_use[1])->getPos();
  q=loop->getAtomAtRes(PID::N,DOF_indices_to_use[2])->getPos();
  b_ang[4]=AngleBetweenVectors(p2-p1,p2-q);

  p1=loop->getAtomAtRes(PID::C,DOF_indices_to_use[1])->getPos();
  p2=loop->getAtomAtRes(PID::N,DOF_indices_to_use[2])->getPos();
  q=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[2])->getPos();
  b_ang[5]=AngleBetweenVectors(p2-p1,p2-q);
  

  p1=loop->getAtomAtRes(PID::N,DOF_indices_to_use[2])->getPos();
  p2=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[2])->getPos();
  q=loop->getAtomAtRes(PID::C,DOF_indices_to_use[2])->getPos();
  b_ang[6]=AngleBetweenVectors(p2-p1,p2-q);
  
  //peptide torsion angles
  //t_ang(1:2) = pi
  //t_ang[0] = pi;
  //t_ang[1] = pi;
  Vector3 q1,q2;
  
  p1=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[0])->getPos();
  p2=loop->getAtomAtRes(PID::C,DOF_indices_to_use[0])->getPos();
  q1=loop->getAtomAtRes(PID::N,DOF_indices_to_use[1])->getPos();
  q2=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[1])->getPos();
  t_ang[0] = PMath::TorsionAngle(p2-p1,q2-q1,q1-p2)*deg2rad;
  if (t_ang[0] < 0)t_ang[0] = t_ang[0] + PI;
  else t_ang[0] = t_ang[0] - PI;

  p1=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[1])->getPos();
  p2=loop->getAtomAtRes(PID::C,DOF_indices_to_use[1])->getPos();
  q1=loop->getAtomAtRes(PID::N,DOF_indices_to_use[2])->getPos();
  q2=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[2])->getPos();
  t_ang[1] = PMath::TorsionAngle(p2-p1,q2-q1,q1-p2)*deg2rad;
  if (t_ang[1] < 0)t_ang[1] = t_ang[1] + PI;
  else t_ang[1] = t_ang[1] - PI;

  initialize_loop_closure(b_len, b_ang, t_ang);
  
  r_n[1][0]=loop->getAtomAtRes(PID::N,DOF_indices_to_use[0])->getPos().x;
  r_n[1][1]=loop->getAtomAtRes(PID::N,DOF_indices_to_use[0])->getPos().y;
  r_n[1][2]=loop->getAtomAtRes(PID::N,DOF_indices_to_use[0])->getPos().z;
  
  r_a[1][0]=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[0])->getPos().x;
  r_a[1][1]=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[0])->getPos().y;
  r_a[1][2]=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[0])->getPos().z;
  
  r_a[3][0] = endPriorG->x;
  r_a[3][1] = endPriorG->y;
  r_a[3][2] = endPriorG->z;
  
  r_c[3][0] = endG->x;
  r_c[3][1] = endG->y;
  r_c[3][2] = endG->z;
  
  //     ! call tripeptide loop closure routine
  //     call solv_3pep_poly(r_n(:,2), r_a(:,2), r_a(:,4), r_c(:,4), &
  //          r_soln_n, r_soln_a, r_soln_c, n_soln)

  int solution_selection = -1;
  solve_3pep_poly(r_n[1], r_a[1], r_a[3], r_c[3], r_soln_n, r_soln_a, r_soln_c, &n_soln, solution_selection);
  Vector3 temp;
  
  double tang;
  Vector3 u0,u1,u2,u3,u4,dest,u,v;
//  IKSolutions solutions;
  solutions.clear();
  if (n_soln >= 1){  
	  //NOTE: Yajia: I changed here. Next one line.
	  assert( solution_selection >= 0);
	  i = solution_selection;
//    for(i=0;i<n_soln;i++)
    {
      IKSolution AllMove;
      ChainMove CMove;
      CMove.blockType = PID::BACKBONE;
      CMove.dir = forward;
      u0 = loop->getAtomAtRes(PID::N,DOF_indices_to_use[0])->getPos();
      u1 = loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[0])->getPos();
      u2 = loop->getAtomAtRes(PID::C,DOF_indices_to_use[0])->getPos();
      u3 = Vector3(r_soln_c[i][0][0],r_soln_c[i][0][1],r_soln_c[i][0][2]);
      tang = PMath::TorsionAngle(u2-u1,u3-u1,u1-u0);
      CMove.DOF_index = DOF_indices_to_use[0]*2;
      CMove.degrees = -tang;
      AllMove.push_back(CMove);
      loop->RotateChain_noGridUpdate(PID::BACKBONE,CMove.DOF_index,forward,CMove.degrees);
      
      u0 = loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[0])->getPos();
      u1 = loop->getAtomAtRes(PID::C,DOF_indices_to_use[0])->getPos();
      u2 = loop->getAtomAtRes(PID::N,DOF_indices_to_use[1])->getPos();
      u3 = Vector3(r_soln_n[i][1][0],r_soln_n[i][1][1],r_soln_n[i][1][2]);
      tang = PMath::TorsionAngle(u2-u1,u3-u1,u1-u0);
      CMove.DOF_index = DOF_indices_to_use[0]*2+1;
      CMove.degrees = -tang;
      AllMove.push_back(CMove);
      loop->RotateChain_noGridUpdate(PID::BACKBONE,CMove.DOF_index,forward,CMove.degrees);
      
      u0 = loop->getAtomAtRes(PID::N,DOF_indices_to_use[1])->getPos();
      u1 = loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[1])->getPos();
      u2 = loop->getAtomAtRes(PID::C,DOF_indices_to_use[1])->getPos();
      u3 = Vector3(r_soln_c[i][1][0],r_soln_c[i][1][1],r_soln_c[i][1][2]);
      tang = PMath::TorsionAngle(u2-u1,u3-u1,u1-u0);
      CMove.DOF_index = DOF_indices_to_use[1]*2;
      CMove.degrees = -tang;
      AllMove.push_back(CMove);
      loop->RotateChain_noGridUpdate(PID::BACKBONE,CMove.DOF_index,forward,CMove.degrees);
      
      u0 = loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[1])->getPos();
      u1 = loop->getAtomAtRes(PID::C,DOF_indices_to_use[1])->getPos();
      u2 = loop->getAtomAtRes(PID::N,DOF_indices_to_use[2])->getPos();
      u3 = Vector3(r_soln_n[i][2][0],r_soln_n[i][2][1],r_soln_n[i][2][2]);
      tang = PMath::TorsionAngle(u2-u1,u3-u1,u1-u0);
      CMove.DOF_index = DOF_indices_to_use[1]*2+1;
      CMove.degrees = -tang;
      AllMove.push_back(CMove);
      loop->RotateChain_noGridUpdate(PID::BACKBONE,CMove.DOF_index,forward,CMove.degrees);
      
      u0 = loop->getAtomAtRes(PID::N,DOF_indices_to_use[2])->getPos();
      u1 = loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[2])->getPos();
      u2 = loop->getAtomAtRes(PID::C,DOF_indices_to_use[2])->getPos();
      u3 = Vector3(r_soln_c[i][2][0],r_soln_c[i][2][1],r_soln_c[i][2][2]);
      tang = PMath::TorsionAngle(u2-u1,u3-u1,u1-u0);
      CMove.DOF_index = DOF_indices_to_use[2]*2;
      CMove.degrees = -tang;
      AllMove.push_back(CMove);
      
      if (endNextG!=NULL){
        loop->RotateChain_noGridUpdate(PID::BACKBONE,CMove.DOF_index,forward,CMove.degrees);
        u0 = loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[2])->getPos();
        u1 = loop->getAtomAtRes(PID::C,DOF_indices_to_use[2])->getPos();
        u2 = loop->getAtomAtRes(PID::O,DOF_indices_to_use[2])->getPos();
        u3 = Vector3(endNextG->x,endNextG->y,endNextG->z);
        tang = PMath::TorsionAngle(u2-u1,u3-u1,u1-u0);
        CMove.DOF_index = DOF_indices_to_use[2]*2+1;
        CMove.degrees = -tang;
        AllMove.push_back(CMove);
        loop->RotateChain_noGridUpdate(PID::BACKBONE,AllMove[4].DOF_index,forward,-AllMove[4].degrees);
      }

      loop->RotateChain_noGridUpdate(PID::BACKBONE,AllMove[3].DOF_index,forward,-AllMove[3].degrees);
      loop->RotateChain_noGridUpdate(PID::BACKBONE,AllMove[2].DOF_index,forward,-AllMove[2].degrees);
      loop->RotateChain_noGridUpdate(PID::BACKBONE,AllMove[1].DOF_index,forward,-AllMove[1].degrees);
      loop->RotateChain_noGridUpdate(PID::BACKBONE,AllMove[0].DOF_index,forward,-AllMove[0].degrees);
      
      solutions.push_back(AllMove);
    }
  }
  return n_soln;
} 

IKSolutions PExactIKSolver::FindSolutions(PProtein *loop, int DOF_indices_to_use[3], Vector3 *endPriorG, Vector3 *endG, Vector3 *endNextG){
//!-----------------------------------------------------------------------
//! This is a sample driver routine to reconstruct tripeptide loops
//! from the coordinates in a pdb file using a canonical bond lengths and
//! angles.
//!-----------------------------------------------------------------------
//  use in_out
//  use tripep_closure
//!-----------------------------------------------------------------------
//  implicit none
//  real(dp) :: time1, time2
//  character(len=100) :: in_pdb, out_prefix, out_pdb
  char in_pdb[100], out_prefix[100], out_pdb[100];
//  character(len=4) :: res_name(5)
  char res_name[5][4];
//  integer :: n0, n_soln, unit = 17, i, j, k, n1, n2, order(max_soln)
  int n0, n_soln, i, j, k, n1, n2;
//  real(dp) :: r_n(3,5), r_a(3,5), r_c(3,5)
  double r_n[5][3], r_a[5][3], r_c[5][3];
//  real(dp) :: r_soln_n(3,3,max_soln), r_soln_a(3,3,max_soln), r_soln_c(3,3,max_soln)
  double r_soln_n[max_soln][3][3], r_soln_a[max_soln][3][3], r_soln_c[max_soln][3][3];
//  real(dp) :: rmsd, sum, dr(3)
  double rmsd, sum, dr[3];
//  real(dp) :: r0_n(3,3), r0_a(3,3), r0_c(3,3)
  double r0_n[3][3], r0_a[3][3], r0_c[3][3];

//  logical :: calc_rmsd, write_out_pdb
  int  calc_rmsd, write_out_pdb;
//  real(dp) :: b_len(6), b_ang(7), t_ang(2)
  double b_len[6], b_ang[7], t_ang[2];

//  ! input parameters (lengths)
//  real(dp), parameter :: b_na = 1.45d0, b_ac = 1.52d0, b_cn = 1.33d0
  double b_na = 1.45e0, b_ac = 1.52e0, b_cn = 1.33e0;
//  ! input parameters (angles in radians)
//  real(dp), parameter :: ang_nac = 111.6*deg2rad, ang_acn = 117.5*deg2rad
  double ang_nac = 112*deg2rad, ang_acn = 120.0*deg2rad;
//  real(dp), parameter :: ang_cna = 120.0*deg2rad
  double ang_cna = 121.9*deg2rad;

//!-----------------------------------------------------------------------

//  call my_timer(time1)

//  !! initialize: set up input parameters
//  ! bond lengths
//  b_len(1:6) = (/ b_ac, b_cn, b_na, b_ac, b_cn, b_na /)
  //cout<<"InFindSolutions2"<<endl;
  Vector3 p1,p2;
  p1=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[0])->getPos();
  p2=loop->getAtomAtRes(PID::C,DOF_indices_to_use[0])->getPos();
  b_len[0]=p1.distance(p2);

  p1=loop->getAtomAtRes(PID::C,DOF_indices_to_use[0])->getPos();
  p2=loop->getAtomAtRes(PID::N,DOF_indices_to_use[1])->getPos();
  b_len[1]=p1.distance(p2);
  p1=loop->getAtomAtRes(PID::N,DOF_indices_to_use[1])->getPos();
  p2=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[1])->getPos();
  b_len[2]=p1.distance(p2);

  p1=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[1])->getPos();
  p2=loop->getAtomAtRes(PID::C,DOF_indices_to_use[1])->getPos();
  b_len[3]=p1.distance(p2);

  p1=loop->getAtomAtRes(PID::C,DOF_indices_to_use[1])->getPos();
  p2=loop->getAtomAtRes(PID::N,DOF_indices_to_use[2])->getPos();
  b_len[4]=p1.distance(p2);

  p1=loop->getAtomAtRes(PID::N,DOF_indices_to_use[2])->getPos();
  p2=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[2])->getPos();
  b_len[5]=p1.distance(p2);


//  ! bond angles
//  b_ang(1:7) = (/ ang_nac, ang_acn, ang_cna, ang_nac, ang_acn, ang_cna, ang_nac /)

  Vector3 q;
  p1=loop->getAtomAtRes(PID::N,DOF_indices_to_use[0])->getPos();
  p2=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[0])->getPos();
  q=loop->getAtomAtRes(PID::C,DOF_indices_to_use[0])->getPos();
  b_ang[0]=AngleBetweenVectors(p2-p1,p2-q);

  p1=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[0])->getPos();
  p2=loop->getAtomAtRes(PID::C,DOF_indices_to_use[0])->getPos();
  q=loop->getAtomAtRes(PID::N,DOF_indices_to_use[1])->getPos();
  b_ang[1]=AngleBetweenVectors(p2-p1,p2-q);

  p1=loop->getAtomAtRes(PID::C,DOF_indices_to_use[0])->getPos();
  p2=loop->getAtomAtRes(PID::N,DOF_indices_to_use[1])->getPos();
  q=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[1])->getPos();
  b_ang[2]=AngleBetweenVectors(p2-p1,p2-q);

  p1=loop->getAtomAtRes(PID::N,DOF_indices_to_use[1])->getPos();
  p2=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[1])->getPos();
  q=loop->getAtomAtRes(PID::C,DOF_indices_to_use[1])->getPos();
  b_ang[3]=AngleBetweenVectors(p2-p1,p2-q);

  p1=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[1])->getPos();
  p2=loop->getAtomAtRes(PID::C,DOF_indices_to_use[1])->getPos();
  q=loop->getAtomAtRes(PID::N,DOF_indices_to_use[2])->getPos();
  b_ang[4]=AngleBetweenVectors(p2-p1,p2-q);

  p1=loop->getAtomAtRes(PID::C,DOF_indices_to_use[1])->getPos();
  p2=loop->getAtomAtRes(PID::N,DOF_indices_to_use[2])->getPos();
  q=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[2])->getPos();
  b_ang[5]=AngleBetweenVectors(p2-p1,p2-q);


  p1=loop->getAtomAtRes(PID::N,DOF_indices_to_use[2])->getPos();
  p2=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[2])->getPos();
  q=loop->getAtomAtRes(PID::C,DOF_indices_to_use[2])->getPos();
  b_ang[6]=AngleBetweenVectors(p2-p1,p2-q);

  //peptide torsion angles
  //t_ang(1:2) = pi
  //t_ang[0] = pi;
  //t_ang[1] = pi;
  Vector3 q1,q2;

  p1=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[0])->getPos();
  p2=loop->getAtomAtRes(PID::C,DOF_indices_to_use[0])->getPos();
  q1=loop->getAtomAtRes(PID::N,DOF_indices_to_use[1])->getPos();
  q2=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[1])->getPos();
  t_ang[0] = PMath::TorsionAngle(p2-p1,q2-q1,q1-p2)*deg2rad;
  if (t_ang[0] < 0)t_ang[0] = t_ang[0] + PI;
  else t_ang[0] = t_ang[0] - PI;

  p1=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[1])->getPos();
  p2=loop->getAtomAtRes(PID::C,DOF_indices_to_use[1])->getPos();
  q1=loop->getAtomAtRes(PID::N,DOF_indices_to_use[2])->getPos();
  q2=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[2])->getPos();
  t_ang[1] = PMath::TorsionAngle(p2-p1,q2-q1,q1-p2)*deg2rad;
  if (t_ang[1] < 0)t_ang[1] = t_ang[1] + PI;
  else t_ang[1] = t_ang[1] - PI;

  initialize_loop_closure(b_len, b_ang, t_ang);

  r_n[1][0]=loop->getAtomAtRes(PID::N,DOF_indices_to_use[0])->getPos().x;
  r_n[1][1]=loop->getAtomAtRes(PID::N,DOF_indices_to_use[0])->getPos().y;
  r_n[1][2]=loop->getAtomAtRes(PID::N,DOF_indices_to_use[0])->getPos().z;

  r_a[1][0]=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[0])->getPos().x;
  r_a[1][1]=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[0])->getPos().y;
  r_a[1][2]=loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[0])->getPos().z;

  r_a[3][0] = endPriorG->x;
  r_a[3][1] = endPriorG->y;
  r_a[3][2] = endPriorG->z;

  r_c[3][0] = endG->x;
  r_c[3][1] = endG->y;
  r_c[3][2] = endG->z;

  //     ! call tripeptide loop closure routine
  //     call solv_3pep_poly(r_n(:,2), r_a(:,2), r_a(:,4), r_c(:,4), &
  //          r_soln_n, r_soln_a, r_soln_c, n_soln)

  int solution_selection = -1;
  solve_3pep_poly(r_n[1], r_a[1], r_a[3], r_c[3], r_soln_n, r_soln_a, r_soln_c, &n_soln, solution_selection);
  Vector3 temp;

  double tang;
  Vector3 u0,u1,u2,u3,u4,dest,u,v;
  IKSolutions solutions;
  solutions.clear();
  if (n_soln >= 1){
	  //NOTE: Yajia: I changed here. Next one line.
	  assert( solution_selection >= 0);
	  i = solution_selection;
//    for(i=0;i<n_soln;i++)
    {
      IKSolution AllMove;
      ChainMove CMove;
      CMove.blockType = PID::BACKBONE;
      CMove.dir = forward;
      u0 = loop->getAtomAtRes(PID::N,DOF_indices_to_use[0])->getPos();
      u1 = loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[0])->getPos();
      u2 = loop->getAtomAtRes(PID::C,DOF_indices_to_use[0])->getPos();
      u3 = Vector3(r_soln_c[i][0][0],r_soln_c[i][0][1],r_soln_c[i][0][2]);
      tang = PMath::TorsionAngle(u2-u1,u3-u1,u1-u0);
      CMove.DOF_index = DOF_indices_to_use[0]*2;
      CMove.degrees = -tang;
      AllMove.push_back(CMove);
      loop->RotateChain_noGridUpdate(PID::BACKBONE,CMove.DOF_index,forward,CMove.degrees);

      u0 = loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[0])->getPos();
      u1 = loop->getAtomAtRes(PID::C,DOF_indices_to_use[0])->getPos();
      u2 = loop->getAtomAtRes(PID::N,DOF_indices_to_use[1])->getPos();
      u3 = Vector3(r_soln_n[i][1][0],r_soln_n[i][1][1],r_soln_n[i][1][2]);
      tang = PMath::TorsionAngle(u2-u1,u3-u1,u1-u0);
      CMove.DOF_index = DOF_indices_to_use[0]*2+1;
      CMove.degrees = -tang;
      AllMove.push_back(CMove);
      loop->RotateChain_noGridUpdate(PID::BACKBONE,CMove.DOF_index,forward,CMove.degrees);

      u0 = loop->getAtomAtRes(PID::N,DOF_indices_to_use[1])->getPos();
      u1 = loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[1])->getPos();
      u2 = loop->getAtomAtRes(PID::C,DOF_indices_to_use[1])->getPos();
      u3 = Vector3(r_soln_c[i][1][0],r_soln_c[i][1][1],r_soln_c[i][1][2]);
      tang = PMath::TorsionAngle(u2-u1,u3-u1,u1-u0);
      CMove.DOF_index = DOF_indices_to_use[1]*2;
      CMove.degrees = -tang;
      AllMove.push_back(CMove);
      loop->RotateChain_noGridUpdate(PID::BACKBONE,CMove.DOF_index,forward,CMove.degrees);

      u0 = loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[1])->getPos();
      u1 = loop->getAtomAtRes(PID::C,DOF_indices_to_use[1])->getPos();
      u2 = loop->getAtomAtRes(PID::N,DOF_indices_to_use[2])->getPos();
      u3 = Vector3(r_soln_n[i][2][0],r_soln_n[i][2][1],r_soln_n[i][2][2]);
      tang = PMath::TorsionAngle(u2-u1,u3-u1,u1-u0);
      CMove.DOF_index = DOF_indices_to_use[1]*2+1;
      CMove.degrees = -tang;
      AllMove.push_back(CMove);
      loop->RotateChain_noGridUpdate(PID::BACKBONE,CMove.DOF_index,forward,CMove.degrees);

      u0 = loop->getAtomAtRes(PID::N,DOF_indices_to_use[2])->getPos();
      u1 = loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[2])->getPos();
      u2 = loop->getAtomAtRes(PID::C,DOF_indices_to_use[2])->getPos();
      u3 = Vector3(r_soln_c[i][2][0],r_soln_c[i][2][1],r_soln_c[i][2][2]);
      tang = PMath::TorsionAngle(u2-u1,u3-u1,u1-u0);
      CMove.DOF_index = DOF_indices_to_use[2]*2;
      CMove.degrees = -tang;
      AllMove.push_back(CMove);

      if (endNextG!=NULL){
        loop->RotateChain_noGridUpdate(PID::BACKBONE,CMove.DOF_index,forward,CMove.degrees);
        u0 = loop->getAtomAtRes(PID::C_ALPHA,DOF_indices_to_use[2])->getPos();
        u1 = loop->getAtomAtRes(PID::C,DOF_indices_to_use[2])->getPos();
        u2 = loop->getAtomAtRes(PID::O,DOF_indices_to_use[2])->getPos();
        u3 = Vector3(endNextG->x,endNextG->y,endNextG->z);
        tang = PMath::TorsionAngle(u2-u1,u3-u1,u1-u0);
        CMove.DOF_index = DOF_indices_to_use[2]*2+1;
        CMove.degrees = -tang;
        AllMove.push_back(CMove);
        loop->RotateChain_noGridUpdate(PID::BACKBONE,AllMove[4].DOF_index,forward,-AllMove[4].degrees);
      }

      loop->RotateChain_noGridUpdate(PID::BACKBONE,AllMove[3].DOF_index,forward,-AllMove[3].degrees);
      loop->RotateChain_noGridUpdate(PID::BACKBONE,AllMove[2].DOF_index,forward,-AllMove[2].degrees);
      loop->RotateChain_noGridUpdate(PID::BACKBONE,AllMove[1].DOF_index,forward,-AllMove[1].degrees);
      loop->RotateChain_noGridUpdate(PID::BACKBONE,AllMove[0].DOF_index,forward,-AllMove[0].degrees);

      solutions.push_back(AllMove);
    }
  }
  return solutions;
}

//returns in radians
double PExactIKSolver::AngleBetweenVectors(Vector3 n1, const Vector3 n2) {
  double d1 = n1.dot(n2)/(n1.norm()*n2.norm());
  if (d1>1) d1 = 1;
  if (d1<-1) d1 = -1;
  return acos(d1);
}
