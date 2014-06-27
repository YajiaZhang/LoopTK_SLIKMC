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

#ifndef PCHAIN_H
#define PCHAIN_H

#include "PLibraries.h"

/*
 * Defines naming conventions and
 * geometric data which apply to
 * all proteins.
 */

#define		NUM_BACKBONE_ATOMS      3		/* N, CA, C are backbone atoms. */
#define		LENGTH_C_N              1.33017		/* C-N bond length in backbone. */
#define		ANGLE_C_N_CA            121.254		/* C-N-CA bond angle in backbone. */
#define		ANGLE_CA_C_N            116.606		/* CA-C-N bond angle in backbone. */

#define		ANGLE_O_C_N	            123.5		/* O=C-N bond angle in backbone. */
#define 	LENGTH_CA_C             1.51
#define		LENGTH_O_C              1.24
#define 	LENGTH_N_CA		1.45

#define   COLLISION_THRESHOLD     0.75
#define PI 3.14159265


/* Macro to check whether a given string identifies a backbone atom. */
#define		isBackbone(x)		((x) == PID::N || (x) == PID::C_ALPHA || (x) == PID::C)
#define		isBackboneOrOxy(x)	(isBackbone(x) || (x) == PID::O)


/**
 * The PID namespace defines standard identifiers for atoms, blocks, residues,
 * and other resources provided with LoopTK.
 */
namespace PID {

  // Atom IDs
  static const string C = "C";
  static const string C_ALPHA = "CA";
  static const string C_BETA = "CB";
  static const string O = "O";
  static const string N = "N";
  static const string S = "S";

  // Block types
  static string BACKBONE = "backbone";
  static string SIDECHAIN = "sidechain";

  // Residues
  static const string ALA = "ALA";
  static const string ARG = "ARG";
  static const string ASN = "ASN";
  static const string ASP = "ASP";
  static const string BACKBONE_RES = "BBR";
  static const string CB_RES = "CBR";
  static const string CYS = "CYS";
  static const string GLN = "GLN";
  static const string GLU = "GLU";
  static const string GLY = "GLY";
  static const string HIS = "HIS";
  static const string ILE = "ILE";
  static const string LEU = "LEU";
  static const string LYS = "LYS";
  static const string MET = "MET";
  static const string PHE = "PHE";
  static const string PRO = "PRO";
  static const string SER = "SER";
  static const string THR = "THR";
  static const string TRP = "TRP";
  static const string TYR = "TYR";
  static const string VAL = "VAL";

  // Residue abbreviations. "_C" stands for "char".
  // These are the standard abbreviations in use now.

  //not used anywhere yet...
  static const char ALA_C = 'A';
  static const char ARG_C = 'R';
  static const char ASN_C = 'N';
  static const char ASP_C = 'D';
  static const char CYS_C = 'C';
  static const char GLN_C = 'Q';
  static const char GLU_C = 'E';
  static const char GLY_C = 'G';
  static const char HIS_C = 'H';
  static const char ILE_C = 'I';
  static const char LEU_C = 'L';
  static const char LYS_C = 'K';
  static const char MET_C = 'M';
  static const char PHE_C = 'F';
  static const char PRO_C = 'P';
  static const char SER_C = 'S';
  static const char THR_C = 'T';
  static const char TRP_C = 'W';
  static const char TYR_C = 'Y';
  static const char VAL_C = 'V';
};

#endif
