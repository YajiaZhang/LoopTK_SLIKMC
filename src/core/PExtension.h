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

/*
 * PExtension defines a variety of auxiliary
 * classes that are not considered part of the
 * "core" LoopTK class hierarchy. These include
 * file I/O classes, initialization classes,
 * test classes, and so forth.
 */

#ifndef PEXTENSION_H
#define PEXTENSION_H

#include "PBasic.h"
#include "PIKAlgorithms.h"
#include "PConstants.h"
// @package Conformation Analysis

#include <libxml/parser.h>
#include <libxml/tree.h>

typedef list<PLightChain *> PLightChainList;

/**
 *
 * Defines a space of protein loop conformations,
 * as well as the static portions of the protein
 * before and after the loop.
 */
class PConformationSpace: public PLightChainList
{
 public:

  /**
   * Constructs a new <code>PConformationSpace</code> with
   * no conformations. The required parameter indicates the
   * static portions of the protein, i.e. subchains before
   * and after the loop; the first <code>PLightChain</code>
   * in the pair comes before the loop, the second comes
   * after it.
   */

  PConformationSpace(pair<PLightChain *,PLightChain *> staticPortions);
  ~PConformationSpace();
  
  /**
   * Gets the two static portions of the loop in the form of PLightChain's.
   * The first PLightChain in the pair comes before the loop; the second
   * comes afterwards.
   */

  pair<PLightChain *, PLightChain *> getStaticPortions() { return m_staticPortions; }

 private:
  pair<PLightChain *,PLightChain *> m_staticPortions;

};


/*
 * A csConfLine (conformation line) defines a new position for one
 * particular atom in one particular loop residue.
 */

#define		csConfLine		pair<pair<int, string>, Vector3>

enum WarningIndicator {
  SUPPRESS_WARNINGS,
  ENABLE_WARNINGS
};
// @package Main Infrastructure
/**
 *
 * Used to initialize the program. A call to <code>Initialize()</code>
 * will load all resources needed by the infrastructure and set
 * warnings on or off.
 */

class LoopTK {
 public:

  /**
   * Initializes LoopTK with the resources found at the directory at path <code>resourcePath</code>.  
   * Sets warnings on or off, as specified by <code>ind</code>.
   */

  static void Initialize(WarningIndicator ind, const string &resourcePath = "resources");

  /**
   * Suppresses all LoopTK warnings.
   */

  static void SuppressWarnings();

  /**
   * Enables all LoopTK warnings.
   */

  static void EnableWarnings();

  /**
   * Displays the specified <code>warning</code>
   * if and only if warnings are enabled.
   */

  static void DisplayWarning(const string &warning);
 private:
  static bool m_warningsEnabled;
};

//@package Resource Management
/**
 *
 * Contains methods to load and initialize resources
 * needed by LoopTK.  This is a private, static class
 * used by <code>LoopTK::Initialize()</code>.
 */

class PInit {
  friend class LoopTK;

  public:
    static string GetName(xmlNodePtr node);
    static string GetAttribute(xmlNodePtr node, const string &attr);

  private:
    static void InitializeResources(const string &resourceDir);

    static HASH_MAP_STRPAIR_EX(vector<PBlockBondShell>) blockBonds;
    static HASH_MAP_STRPAIR_EX(HASH_MAP_STRPAIR_EX(Vector3)) relPos;

    static void setupBlockShells(const string &file);
    static void setupBlockConnections(const string &file);
    static void setupResidueShells(const string &file);

    static void setupAtomShells(const string &file);
    static void setupChiIndices(const string &file);
    static void setupEpsilonVals(const string &file);
    static void setupIDMaps(const string &file);
    static void setupRotamers(const string &file);
};
//@package Input-Output
/**
 * 
 *
 * Contains methods to read PChains from a PDB file
 * and write them back out in PDB format. See the
 * <a href="http://www.pdb.org">PDB homepage</a>
 * for more information.
 */

class PDBIO {
  friend class CS2IO;	/* To access readFromLines(). */

  public:

    /**
     * A convenience structure that encapsulates
     * the data on one line of a PDB "ATOM" descriptor.
     */

    struct pdbAtomLine {
      string atomName, resName, elem;
      char altLoc, insCode, chainID;
      int atomNum, resNum;
      string occupancy, tempFactor;
      Vector3 pos;
    };

    /**
     * Reads the specified chain in PDB file <code>/fileName</code> 
     * and returns the contents as a PProtein pointer, or NULL if an
     * I/O error is encountered. If <code>chainId</code> is unspecified,
     * the first chain (terminated by "TER") will be read.
     */

    static PProtein* readFromFile(const string &fileName, const string &chainId="");

    /**
     * Writes the specified PChain <code>chain</code> to a file named
     * <code>fileName</code>, using the standard PDB format.  Optionally,
     * a <code>prefix</code> may be appended to each line of output.
     */

    static void writeToFile(PChain *chain, const string &fileName, const string &prefix = "");

    /**
     * Returns <code>true</code> if the PDB file specified by
     * <code>fileName</code> appears to represent a cyclic protein,
     * <code>false</code> otherwise.
     */

    static bool isCyclicPDB(const string &fileName);

    /**
     * Returns a <code>vector</code> containing all the lines in
     * <code>pdbLines</code> that begin with <code>"ATOM  "</code>
     * or <code>"HETATM"</code>, up to the <code>"TER"</code>
     * terminating line (if one is present).
     */

    static vector<string> getAllAtomLines(const vector<string> &pdbLines);

    /**
     * Returns a pdbAtomLine structure containing
     * data from the specified line in a PDB file.
     */

    static pdbAtomLine parseAtomLine(const string &line);

  private:

    /* Helper functions for readFromFile(). */

    static PProtein* readFromLines(const vector<string> &pdbLines);
    static PProtein* readFromMap(const hash_map<int, vector<pdbAtomLine> > &atomMap, 
				 const set<int> &allResNums);
    static void updateProtein(PProtein* &protein, const string &resName,
				PResidueSpec &spec, int resNum);

    static void trimHeadAndTail(vector<int> &resIndices, const hash_map<int,
				vector<pdbAtomLine> > &atomMap);
    static PResidueSpec getSpec(const vector<pdbAtomLine> &atomLines);
    static bool definesBackbone(PResidueSpec &spec);

    /* Helper functions for writeToFile() that do some error checking. */

    static void addAtomNum(string &pdbOutLine, int atomNum, string::size_type prefixLen);
    static void addAtomName(string &pdbOutLine, const string &atomName, string::size_type prefixLen);
    static void addResName(string &pdbOutLine, const string &resName, string::size_type prefixLen);
    static void addChainID(string &pdbOutLine, char chainID, string::size_type prefixLen);
    static void addResNum(string &pdbOutLine, int resNum, string::size_type prefixLen);
    static void addInsertionCode(string &pdbOutLine, const string &insCode, string::size_type prefixLen);
    static void addAtomPos(string &pdbOutLine, const Vector3 &atomPos, string::size_type prefixLen);
    static void addAuxData(string &pdbOutLine, Real occupancy, Real tempFactor,
				const string &segmentID, string::size_type prefixLen);
    static void addElemName(string &pdbOutLine, const string &elemName, string::size_type prefixLen);

};

//@package Math
/**
 * 
 *
 * Contains mathematical functions needed in various
 * LoopTK computations, such as cyclic coordinate descent
 * (CCD) and chain manipulations.
 */

class PMath {
  public:

   /**
    * Computes the angle between the planes determined by the specified four
    * points.
    */

   static Real AngleBetweenPlanes(const Vector3 &p1, const Vector3 &p2,
                                  const Vector3 &p3, const Vector3 &p4);

   /**
    * Computes the angle (in radians) between vectors
    * v1 and v2.
    */

    static Real getAngle(const Vector3 &v1, const Vector3 &v2);

   /**
    * Computes the position of the point "next" 
    * given points prev and cur, the angle prev-cur-next,
    * and the distance between cur and next.
    */

    static Vector3 ComputePos(const Vector3 &prev, const Vector3 &cur, Real angle, Real bondLength);

    /**
     * Returns <code>true</code> if the sphere with specified
     * center and radius intersects any part of the cube with
     * the specified pair of opposing vertices.
     */

    static bool sphereIntersectsCube(const Vector3 &sphereCenter, Real sphereRadius,
				     const Vector3 &cubeVertex1, const Vector3 &cubeVertex2);

    /**
     * Given the vector v, computes and returns
     * a vector not parallel to v.
     */

    static Vector3 ArbitraryNonParallel(const Vector3 &v);

    /**
     * Divides each element in v by the sum of all
     * elements in v.
     */
    static void normalize(vector<double> &v);

    //finds transform that moves s1 to e1 and s2 to e2
    static Matrix4 Find2PtTransform(Vector3 s1, Vector3 s2, Vector3 e1, Vector3 e2);

    /**
     * Finds and returns a rotation matrix for
     * rotating a vector the specified number
     * of radians around the specified axis according to left-hand rule.
     */

    static Matrix3 FindRotationMatrix(const Vector3 &axis, Real angle);

    /**
     * Converts the vector of defined relative
     * coordinates into a vector of global coordinates.
     */

    static Vector3 LocalToGlobalCoord(const vector<Vector3> definedRelativeCoordinates, const vector<Vector3> definedGlobalCoordinates);
    /**
    Converts a 3X3 rotation matrix specified by <code> input </code> to a quaternion vector returned in <code>output</code>.
    */
    static void ConvertToQuaternion(double** input, Vector4* output, vector<int> sign);
    /**
    Interpolates between two quaternions specified by <code>qi</code> and <code>qf</code> and returns the interpolated quternion in <code>qu</code>. <code>u</code> ranges from 0.0 to 1.0
    */
    static void InterpolateQuaternion(Vector4 qi, Vector4 qf, Vector4* qu, double u);
    /**
    Converts a change in quaternion space specified by difference between <code>qi</code> and <code>qf</code> to a change in velocity space and returns it in <code>deltaV</code>.
    */
    static void DeltaQuaterToDeltaVel(Vector4 qi, Vector4 qf, Vector3* delatV);
   
   /** Returns the signum of <code>a</code>*/
     static int signum(double a);
  
    /**returns the torsion angle between the plane defined by vectors <code>w</code> and <code>u</code> and the plane defined by vectors <code>w</code> and <code>v</code>.*/
    static double TorsionAngle(Vector3 u,Vector3 v,Vector3 w);

    static double abs(double d);
};

/**
 * The <code>Angle</code> class is designed to simplify, reduce
 * confusion, and provide consistency in angle-related
 * computations. It eliminates the need to remember
 * whether a function parameter or variable is expressed
 * in radians or degrees.
 */

class Angle {

  public:

    /**
     * Constructs a new <code>Angle</code> of the specified
     * <code>type</code> (<code>angleDeg</code> or <code>angleRad</code>)
     * and with the specified <code>value</code>.
     */

    Angle(AngleType type, Real value);

    /**
     * Returns the value of this <code>Angle</code> in degrees.
     */

    Real getDegrees();

    /**
     * Returns the value of this <code>Angle</code> in radians.
     */

    Real getRadians();

  private:

    AngleType m_type;	/* The angle type (radians or degrees). */
    Real m_value;	/* The angle value in m_type units.	*/

};
//@package Input-Output
/**
 * 
 *
 * Contains methods to read and write
 * from CS2 (second-generation chainspace)
 * files. See DATA_COLLECTOR_DESIGN for
 * more information about these files.
 */

class CS2IO {

  public:

    /**
     * Reads the <code>confNum</code>'th loop conformation from
     * <code>cs2File</code>.  This function creates a finalized
     * version of the entire protein: the <code>PProtein</code>
     * pointer returned to the caller points to the loop, and
     * the entire chain may be accessed with the method
     * <code>getTopLevelProtein()</code>. Conformation "1" is the 
     * original conformation from when this file was created.
     */

    static PProtein* readConformation(const string &cs2File, int confNum);

    /**
     * Reads the <code>confNum</code>'th loop conformation from
     * <code>cs2File</code>.  This function returns a
     * <code>PLightChain</code> containing only data for
     * the loop, with no information about the protein's static
     * portions.
     */

    static PLightChain *readConformationLoopOnly(const string &cs2File, int confNum);
    
    /**
     * Reads the entire conformation space contained in
     * <code>cs2File</code>, including the protein's static
     * portions and all loop conformations.  This function
     * allocates dynamic memory, so the resulting
     * <code>PConformationSpace</code> should be deleted by
     * the caller when no longer needed.
     */

    static PConformationSpace *readConformationSpace(const string &cs2File);

    /**
     * Reads only the static portions of the protein
     * represented in <code>cs2File</code>.  The first
     * member of the <code>pair</code> returned is the
     * subchain before the loop; the second member is
     * the subchain after the loop.
     */

    static pair<PLightChain *, PLightChain *> getNonLoopChains(const string &cs2File);

    /**
     * Appends a new conformation to the specified CS2 file.
     * <code>loop</code> is assumed to represent only that subchain
     * of the parent chain whose conformation has changed.  Furthermore,
     * the file <code>cs2File</code> must already contain at least one
     * conformation of <code>loop</code> or the program will abort.
     */
   
    static void appendConformation(const string &cs2File, PProtein *loop);

    /**
     * Writes a new CS2 file to the file specified by <code>cs2file</code>.
     * The output is a new file containing the PDB representation of <code>loop</code>'s
     * parent protein; the single "loop" line indicating the start
     * and end residues represented by <code>loop</code>; and the single
     * conformation of the atoms in <code>loop</code>.
     */

    static void writeToFile(const string &cs2File, PProtein *loop);

    /**
     * Returns the number of conformations present within the given file.
     */

    static int getNumConformations(const string &cs2File);

    /**
     * Creates a directory of PDB files given input <code>cs2File</code>,
     * one for each conformation.
     */

    static void createPDBDirectory(const string &cs2File,const string &directoryName);


  private:

    /* Versions of the public functions that read from constituent components
     * of CS2 file. */
    static PLightChain *readConformationLoopOnly(const vector<string> &atomLines,
                                                 const vector<string> &loopLines,
                                                 const vector<string> &confLines,
                                                 const string &fileName);
    static pair<PLightChain *, PLightChain *> getNonLoopChains(const vector<string> &atomLines,
                                                               const vector<string> &loopLines,
                                                                const string &fileName);

    /* Grabs the PDB atom lines, "loop" line and (optionally) conformation lines in cs2File. */
    static void parseFile(const string &cs2File, vector<string> &atomLines, vector<string> &loopLines);
    static void parseFile(const string &cs2File, vector<string> &atomLines, vector<string> &loopLines,
          vector<vector<string> > &allConfLines);
    static void parseFile(const string &cs2File, int confNum, vector<string> &atomLines, vector<string> &loopLines,
          vector<string> &confLines);

    /* Helpers for parseFile. */
    static vector<string> extractAtomLines(const vector<string> &cs2Lines, const string &fileName);
    static vector<string> extractLoopLines(const vector<string> &cs2Lines, const string &fileName);
    static vector<vector<string> > extractConfLines(const vector<string> &cs2Lines, const string &fileName);

    /* Parses a "cn" line, using sanityCheck to make sure residue indices are valid. */
    static csConfLine parseConfLine(const string &cs2File, const string &cs2Line,
        const pair<int, int> &sanityCheck);

    /* Counts the number of conformations in the specified CS2 file. */
    static int countConformations(const string &cs2File);

    /* Writes loop to the specified file as the confNum'th conformation,
     * using the specified start and end residue indices into the top level. */
    static void outputHelper(const string &cs2File, PProtein *loop, int confNum,
           int startIndex, int endIndex);

    /* Parses a "loop" line from the specified file into its constituent pair of ints. */
    static pair<int, int> parseLoopLine(const string &cs2File, const string &loopLine);

    /* Applies the specified conformation to loop, subtracting offset from start index. */
    static void applyConfLine(PProtein* &loop, const csConfLine &conf, int offset = 0);
};

typedef Real(*EnergyCalcFn)(const PAtom *, const PAtom *);

//@package Math
/**
 * 
 *
 * <code>PEnergy</code> contains methods to calculate
 * the free energy for pairs of atoms and proteins.
 */

class PEnergy {

  public:

	/**
	 * The default distance threshold used by <code>energyOfChain</code>
	 * if none is provided.
	 */

	static const Real DEFAULT_THRESHOLD = 10;

	/**
	 * Calculates the van der Waals energy between atoms
	 * <code>a1</code> and <code>a2</code>.
	 */

	static Real vanDerWaalsEnergy(const PAtom *a1, const PAtom *a2);

	/**
	 * Calculates the potential energy for the entire chain
	 * specified by <code>chain</code>, ignoring pairs of
         * bonded atoms as well as pairs of atoms more than
         * <code>threshold</code> Angstroms apart.  A user-defined
         * <code>EnergyCalcFn</code> may be specified to define
	 * the potential energy of a pair of atoms, or the default
	 * <code>vanDerWaalsEnergy</code> will be used.
	 */

	static Real energyOfChain(PChain *chain, EnergyCalcFn energyFn = vanDerWaalsEnergy, Real threshold = DEFAULT_THRESHOLD);

	/**
	 * Calculates the collision energy between atoms 
	 * <code>a1</code> and <code>a2</code>. If the radius of 
	 * <code>a1</code> and <code>a2</code> are r1 and r2 respectively,
	 * the collision energy is 1/d^2 - 1/d0^2, where d is the distance
	 * between two atom centers, and d0=COLLISION_THRESHOLD*(r1+r2).
	 * The COLLISION_THRESHOLD is set to 0.75 in PConstants.h.
	 */

	static Real collisionEnergy(const PAtom *a1, const PAtom *a2);


  private:

};
 //@package Optimization
/**
 *
 *
 * Contains methods to Manipulate atoms of multiple loops of
 * same protein or multiple proteins 
 */
class PMoveAtom {

  public:
   /**
   * Returns vector of <code>IKSolutions</code> corresponding to the loops after each 
   * atoms at positions 
   * specified by vector <code>AtomToMove</code> undergo displacements specified by
   * <code>MovementDirn</code>. <code>Dofs</code> are specified for each loop
   * <code>IKSolutions</code> are applied to the <code>loops</code>.
   */

	static vector<IKSolutions> MoveAtom(vector<PProtein*> loops, vector<vector<CDof> > Dofs, vector<int> AtomToMove, vector<Vector3> MovementDirn);
   /**
   * Returns <code>IKSolutions</code> corresponding to the loop after atom 
   * at position <code>AtomToMove</code> undergoes a displacement specified by
   * <code>MovementDirn</code>. All backbone DOFs are used.
   * <code>IKSolutions</code> is applied to the <code>loop</code>.
   */
	static IKSolutions MoveAtom(PProtein *loop, int AtomToMove, Vector3 MovementDirn);

  private:

};

#endif
