/* ===============================================================================================
   AlphaMol: ia program for computing geometric measures of a union of balls

   Author:  Patrice Koehl (collaboration with Herbert Edelsbrunner)
   Date:    9/22/2019
   Version: 1

   Copyright (C) 2019 Patrice Koehl

   This program and all its dependencies are free software; 
   you can redistribute it and/or modify it under the terms 
   of the GNU Lesser General Public License as published by 
   the Free Software Foundation; either version 2.1 of the License, 
   or (at your option) any later version.

   This software is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this software; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

   =============================================================================================== */

/* ===============================================================================================
   Static variables and Local Includes
   =============================================================================================== */

#include <sys/time.h>
#include "AlphaMol.h"
#include "CheckDeriv.h"

/* ===============================================================================================
   Main program
   =============================================================================================== */

int main(int argc, char **argv)
{

/*	==========================================================================================
	Show usage if needed
	========================================================================================== */

	if( argc < 2 )
	{
		usage(argv);
		return -1;
	}

	std::string input = argv[1];
	if( input == "-h" || input == "-help" )
	{
		usage(argv);
		return -1;
	}

/*	==========================================================================================
	Read in all inputs (some values may have been preset)
	========================================================================================== */

	std::string INfile;
	std::string OUTfile;
	int flag_CA;
	double r_h2o = 1.4;
	int flag_deriv = 0;
    double delaunay_eps = 1;

        if (!parse_args(argc, argv, &INfile, &flag_CA, &r_h2o, &flag_deriv,
		&OUTfile)) return 1;

/*	==========================================================================================
	Read in molecule
	========================================================================================== */
	ReadInput read;

	std::vector<Atoms> atoms;

	std::size_t found = INfile.find("crd");

	if(found !=std::string::npos) {
		read.readFromCRD(INfile, atoms, r_h2o);
	} else {
		found = INfile.find("pqr");
		if(found !=std::string::npos) {
			read.readFromPQR(INfile, flag_CA, atoms, r_h2o);
		} else {
			found = INfile.find("pdb");
			if(found !=std::string::npos) {
				read.readFromPDB(INfile, flag_CA, atoms, r_h2o);
			} else {
				std::cout << " " << std::endl;
				std::cout << "Input file format not recognized; program can only read CRD, PQR, and PDB files" << std::endl;
				std::cout << " " << std::endl;
				exit(1);
			}
		}
	}

	std::cout << " " << std::endl;
	std::cout << "Input file                : " << INfile << std::endl;
	std::cout << "Number of atoms (balls)   : " << atoms.size() << std::endl;
	std::cout << "Probe radius              : " << r_h2o << std::endl;
	std::cout << " " << std::endl;

/*	==========================================================================================
	Compute Delaunay triangulation
	========================================================================================== */

	clock_t start_s, stop_s;

	std::vector<Vertex> vertices;
	std::vector<Tetrahedron> tetra;

	int natoms = atoms.size();
	double *coord = new double[3*natoms];
	double *radii = new double[natoms];
	double *coefS = new double[natoms];
	double *coefV = new double[natoms];
	double *coefM = new double[natoms];
	double *coefG = new double[natoms];

	for(int i = 0; i < natoms; i++) {
		for(int j = 0; j < 3; j++) coord[3*i+j] = atoms[i].Coordinates[j];
		radii[i] = atoms[i].Radius;

		coefS[i] = 1.0;
		coefV[i] = 1.0;
		coefM[i] = 1.0;
		coefG[i] = 1.0;
	}
	
	delcx.setup(natoms, coord, radii, coefS, coefV, coefM, coefG, vertices, tetra);

	start_s = clock();
	delcx.regular3D(vertices, tetra, delaunay_eps);
	stop_s = clock();
	std::cout << "Delaunay compute time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " seconds" << std::endl;

/*	==========================================================================================
	Generate alpha complex (with alpha=0.0)
	========================================================================================== */

	start_s = clock();
	double alpha = 0;
	alfcx.alfcx(alpha, vertices, tetra);
	stop_s = clock();
	std::cout << "AlphaCx compute time : " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " seconds" << std::endl;

/*	==========================================================================================
	Compute surface area and, optionally volume of the union of balls.
	If requested, compute also their derivatives
	========================================================================================== */

	std::vector<Edge> edges;
	std::vector<Face> faces;
	alfcx.alphacxEdges(tetra, edges);
	alfcx.alphacxFaces(tetra, faces);

	double Surf, WSurf, Vol, WVol, Mean, WMean, Gauss, WGauss;

	int nfudge = 8;
	double *ballwsurf = new double[natoms+nfudge];
	double *dsurf = new double[3*(natoms+nfudge)];
	memset(dsurf, 0, 3*(natoms+nfudge)*sizeof(double));

	double *ballwvol, *dvol;
	ballwvol = new double[natoms+nfudge];
	dvol = new double[3*(natoms+nfudge)];
	memset(dvol, 0, 3*(natoms+nfudge)*sizeof(double));

	double *ballwmean, *dmean;
	ballwmean = new double[natoms+nfudge];
	dmean = new double[3*(natoms+nfudge)];
	memset(dmean, 0, 3*(natoms+nfudge)*sizeof(double));

	double *ballwgauss, *dgauss;
	ballwgauss = new double[natoms+nfudge];
	dgauss = new double[3*(natoms+nfudge)];
	memset(dgauss, 0, 3*(natoms+nfudge)*sizeof(double));

	start_s = clock();
	volumes.ball_dvolumes(vertices, tetra, edges, faces, &WSurf, &WVol,
	&WMean, &WGauss, &Surf, &Vol, &Mean, &Gauss, ballwsurf, ballwvol,
	ballwmean, ballwgauss, dsurf, dvol, dmean, dgauss, flag_deriv);
	stop_s = clock();

	std::cout << "Volumes compute time : " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " seconds" << std::endl;

	std::cout << " " << std::endl;
	std::cout << "Biomolecule from file      : " << INfile << std::endl;
	std::cout << "Number of atoms (balls)    : " << natoms << std::endl;
	std::cout << "Probe radius               : " << r_h2o << std::endl;
	std::cout << "Unweighted surface area    : " << std::setw(16) << std::fixed << std::setprecision(8) << Surf << std::endl;
	std::cout << "Weighted surface area      : " << std::setw(16) << std::fixed << std::setprecision(8) << WSurf << std::endl;
	std::cout << "Unweighted volume          : " << std::setw(16) << std::fixed << std::setprecision(8) << Vol << std::endl;
	std::cout << "Weighted volume            : " << std::setw(16) << std::fixed << std::setprecision(8) << WVol << std::endl;
	std::cout << "Unweighted mean curvature  : " << std::setw(16) << std::fixed << std::setprecision(8) << Mean << std::endl;
	std::cout << "Weighted mean curvature    : " << std::setw(16) << std::fixed << std::setprecision(8) << WMean << std::endl;
	std::cout << "Unweighted Gauss curvature : " << std::setw(16) << std::fixed << std::setprecision(8) << Gauss << std::endl;
	std::cout << "Weighted Gauss curvature   : " << std::setw(16) << std::fixed << std::setprecision(8) << WGauss << std::endl;
	std::cout << " " << std::endl;

	if(flag_deriv == 1) {

		std::cout << "Compare analytical with numerical derivatives: " << std::endl;
		CheckDeriv(natoms, coord, radii, coefS, coefV, coefM, coefG,
		dsurf, dvol, dmean, dgauss);

	}

	delete [] coord; delete [] radii; delete [] coefS; delete [] coefV; delete [] coefM; delete [] coefG;
	delete [] ballwsurf; delete [] dsurf; delete [] ballwvol; delete [] dvol;
	delete [] ballwmean; delete [] dmean; delete [] ballwgauss; delete [] dgauss;

	return 0;

}

/* ===============================================================================================
   Usage
   =============================================================================================== */

static void usage(char** argv)
{
    std::cout << "\n\n" <<std::endl;
    std::cout << "     " << "================================================================================================"<<std::endl;
    std::cout << "     " << "================================================================================================"<<std::endl;
    std::cout << "     " << "=                                                                                              ="<<std::endl;
    std::cout << "     " << "=                                         AlphaMol                                             ="<<std::endl;
    std::cout << "     " << "=                                                                                              ="<<std::endl;
    std::cout << "     " << "=     This program computes geometric properties of a molecule represented by a union of balls ="<<std::endl;
    std::cout << "     " << "=                                                                                              ="<<std::endl;
    std::cout << "     " << "=     Usage is:                                                                                ="<<std::endl;
    std::cout << "     " << "=          AlphaMol -i INFILE -r r_h2o - c flag_ca -o OUTFILE                                  ="<<std::endl;
    std::cout << "     " << "=     where:                                                                                   ="<<std::endl;
    std::cout << "     " << "=                 -i INFILE       --> Input file (CRD, PDB, or PQR file)                       ="<<std::endl;
    std::cout << "     " << "=                 -o OUTFILE      --> Output file                                              ="<<std::endl;
    std::cout << "     " << "=                 -r r_h2o        --> probe radius (default 1.4)                               ="<<std::endl;
    std::cout << "     " << "=                 -d deriv        --> flag: no derivatives (0), or derivatives (1)             ="<<std::endl;
    std::cout << "     " << "=                 -c flag_ca      --> For proteins, CA only (0) or all atoms (1) (default 1)   ="<<std::endl;
    std::cout << "     " << "================================================================================================"<<std::endl;
    std::cout << "     " << "================================================================================================"<<std::endl;
    std::cout << "\n\n" <<std::endl;
}

/* ===============================================================================================
   Parse Argument from command line:

   =============================================================================================== */

bool parse_args(int argc, char **argv, std::string *INfile, int *flag_ca, double *r_h2o, 
	int *flag_deriv, std::string *OUTfile) 
{

//
// Make sure we have at least two parameters....
//
	std::string param;
	if (argc == 1)
	{
		return false;
	}
	else
	{
		for (int i = 1; i < argc - 1; i = i + 2)
		{
			param = argv[i];

			if (param == "-i") {
				*INfile = argv[i + 1];
			}
			else if (param == "-o") {
				*OUTfile = argv[i + 1];
			}
			else if(param == "-c") {
				*flag_ca = std::atoi(argv[i + 1]);
			}
			else if(param == "-r") {
				*r_h2o = std::atof(argv[i + 1]);
			}
			else if(param == "-d") {
				*flag_deriv = std::atoi(argv[i + 1]);
			}
		}
  	}
	return true;
}



JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    mod.method("get_geometric_measures",  [](
        jlcxx::ArrayRef<double> outs, 
        const jlcxx::ArrayRef<double> in_coordinates, 
        const jlcxx::ArrayRef<double> in_radii, 
        const double probe_radius,
        const double delaunay_eps){

        /*	==========================================================================================
            Compute Delaunay triangulation
            ========================================================================================== */

        clock_t start_s, stop_s;

        std::vector<Vertex> vertices;
        std::vector<Tetrahedron> tetra;

        int natoms = in_coordinates.size()/3;
        double *coord = new double[3*natoms];
        double *radii = new double[natoms];
        double *coefS = new double[natoms];
        double *coefV = new double[natoms];
        double *coefM = new double[natoms];
        double *coefG = new double[natoms];



        for(int i = 0; i < natoms; i++) {
            for(int j = 0; j < 3; j++){
                coord[3*i+j] = in_coordinates[3*i+j];
            }
            radii[i] = in_radii[i] + probe_radius;
            coefV[i] = 1.0;
            coefS[i] = 1.0;
            coefM[i] = 1.0;
            coefG[i] = 1.0;
        }

        delcx.setup(natoms, coord, radii, coefS, coefV, coefM, coefG, vertices, tetra);

        delcx.regular3D(vertices, tetra, delaunay_eps);

        /*	==========================================================================================
            Generate alpha complex (with alpha=0.0)
            ========================================================================================== */

        double alpha = 0;
        alfcx.alfcx(alpha, vertices, tetra);

        /*	==========================================================================================
            Compute surface area and, optionally volume of the union of balls.
            If requested, compute also their derivatives
            ========================================================================================== */

        std::vector<Edge> edges;
        std::vector<Face> faces;
        alfcx.alphacxEdges(tetra, edges);
        alfcx.alphacxFaces(tetra, faces);

        double Surf, WSurf, Vol, WVol, Mean, WMean, Gauss, WGauss;

        int nfudge = 8;
        double *ballwsurf = new double[natoms+nfudge];
        double *dsurf = new double[3*(natoms+nfudge)];
        memset(dsurf, 0, 3*(natoms+nfudge)*sizeof(double));

        double *ballwvol, *dvol;
        ballwvol = new double[natoms+nfudge];
        dvol = new double[3*(natoms+nfudge)];
        memset(dvol, 0, 3*(natoms+nfudge)*sizeof(double));

        double *ballwmean, *dmean;
        ballwmean = new double[natoms+nfudge];
        dmean = new double[3*(natoms+nfudge)];
        memset(dmean, 0, 3*(natoms+nfudge)*sizeof(double));

        double *ballwgauss, *dgauss;
        ballwgauss = new double[natoms+nfudge];
        dgauss = new double[3*(natoms+nfudge)];
        memset(dgauss, 0, 3*(natoms+nfudge)*sizeof(double));

        volumes.ball_dvolumes(vertices, tetra, edges, faces, &WSurf, &WVol,
                              &WMean, &WGauss, &Surf, &Vol, &Mean, &Gauss, ballwsurf, ballwvol,
                              ballwmean, ballwgauss, dsurf, dvol, dmean, dgauss, 0);

        outs[0] = Vol;
        outs[1] = Surf;
        outs[2] = Mean;
        outs[3] = Gauss;

        delete [] coord; delete [] radii; delete [] coefS; delete [] coefV; delete [] coefM; delete [] coefG;
        delete [] ballwsurf; delete [] dsurf; delete [] ballwvol; delete [] dvol;
        delete [] ballwmean; delete [] dmean; delete [] ballwgauss; delete [] dgauss;
    });

    mod.method("get_geometric_measures_and_overlap_value",  [](
        jlcxx::ArrayRef<double> outs, 
        const jlcxx::ArrayRef<double> in_coordinates,
        const int molecule_size, 
        const jlcxx::ArrayRef<double> in_radii, 
        const double probe_radius,
        const double overlap_existence_penalty,
        const double overlap_penalty_slope,
        const double delaunay_eps
        ){

        /*	==========================================================================================
            Compute Delaunay triangulation
            ========================================================================================== */

        clock_t start_s, stop_s;

        std::vector<Vertex> vertices;
        std::vector<Tetrahedron> tetra;

        int natoms = in_coordinates.size()/3;
        double *coord = new double[3*natoms];
        double *radii = new double[natoms];
        double *coefS = new double[natoms];
        double *coefV = new double[natoms];
        double *coefM = new double[natoms];
        double *coefG = new double[natoms];


        for(int i = 0; i < natoms; i++) {
            for(int j = 0; j < 3; j++){
                coord[3*i+j] = in_coordinates[i*3+j];
            }
            radii[i] = in_radii[i] + probe_radius;
            coefV[i] = 1.0;
            coefS[i] = 1.0;
            coefM[i] = 1.0;
            coefG[i] = 1.0;
        }

        delcx.setup(natoms, coord, radii, coefS, coefV, coefM, coefG, vertices, tetra);

        delcx.regular3D(vertices, tetra, delaunay_eps);

        /*	==========================================================================================
            Generate alpha complex (with alpha=0.0)
            ========================================================================================== */

        double alpha = 0;
        alfcx.alfcx(alpha, vertices, tetra);

        /*	==========================================================================================
            Compute surface area and, optionally volume of the union of balls.
            If requested, compute also their derivatives
            ========================================================================================== */

        std::vector<Edge> edges;
        std::vector<Face> faces;
        alfcx.alphacxEdges(tetra, edges);
        alfcx.alphacxFaces(tetra, faces);

        double Surf, WSurf, Vol, WVol, Mean, WMean, Gauss, WGauss;

        int nfudge = 8;
        double *ballwsurf = new double[natoms+nfudge];
        double *dsurf = new double[3*(natoms+nfudge)];
        memset(dsurf, 0, 3*(natoms+nfudge)*sizeof(double));

        double *ballwvol, *dvol;
        ballwvol = new double[natoms+nfudge];
        dvol = new double[3*(natoms+nfudge)];
        memset(dvol, 0, 3*(natoms+nfudge)*sizeof(double));

        double *ballwmean, *dmean;
        ballwmean = new double[natoms+nfudge];
        dmean = new double[3*(natoms+nfudge)];
        memset(dmean, 0, 3*(natoms+nfudge)*sizeof(double));

        double *ballwgauss, *dgauss;
        ballwgauss = new double[natoms+nfudge];
        dgauss = new double[3*(natoms+nfudge)];
        memset(dgauss, 0, 3*(natoms+nfudge)*sizeof(double));

        volumes.ball_dvolumes(vertices, tetra, edges, faces, &WSurf, &WVol,
                              &WMean, &WGauss, &Surf, &Vol, &Mean, &Gauss, ballwsurf, ballwvol,
                              ballwmean, ballwgauss, dsurf, dvol, dmean, dgauss, 0);

        double total_overlap_penalty = 0.0;
        
        for (int i = 0; i < edges.size(); i++) {
            int ia = edges[i].Vertices[0] - 4;
            int ib = edges[i].Vertices[1] - 4;

            // Check if the edge is between two different molecules
            if (ia / molecule_size != ib / molecule_size) {
                double edge_ol = in_radii[ia] + in_radii[ib] - edges[i].Length;
                if (edge_ol > 0.0) {
                    total_overlap_penalty += overlap_existence_penalty + (overlap_penalty_slope * edge_ol);
                }
            }
        }

        outs[0] = Vol;
        outs[1] = Surf;
        outs[2] = Mean;
        outs[3] = Gauss;
        outs[4] = total_overlap_penalty;

        delete [] coord; delete [] radii; delete [] coefS; delete [] coefV; delete [] coefM; delete [] coefG;
        delete [] ballwsurf; delete [] dsurf; delete [] ballwvol; delete [] dvol;
        delete [] ballwmean; delete [] dmean; delete [] ballwgauss; delete [] dgauss;
    });

    mod.method("get_geometric_measures_with_derivatives",  [](
            jlcxx::ArrayRef<double> measures_out,
            jlcxx::ArrayRef<double> dvol_out,
            jlcxx::ArrayRef<double> dsurf_out,
            jlcxx::ArrayRef<double> dmean_out,
            jlcxx::ArrayRef<double> dgauss_out,
            const jlcxx::ArrayRef<double> in_coordinates, 
            const jlcxx::ArrayRef<double> in_radii, 
            const double probe_radius,
            const double delaunay_eps){


        /*	==========================================================================================
            Compute Delaunay triangulation
            ========================================================================================== */

        clock_t start_s, stop_s;

        std::vector<Vertex> vertices;
        std::vector<Tetrahedron> tetra;

        int natoms = in_coordinates.size()/3;

        double *coord = new double[3*natoms];
        double *radii = new double[natoms];
        double *coefS = new double[natoms];
        double *coefV = new double[natoms];
        double *coefM = new double[natoms];
        double *coefG = new double[natoms];



        for(int i = 0; i < natoms; i++) {
            for(int j = 0; j < 3; j++){
                coord[3*i+j] = in_coordinates[3*i+j];
            }
            radii[i] = in_radii[i] + probe_radius;
            coefV[i] = 1.0;
            coefS[i] = 1.0;
            coefM[i] = 1.0;
            coefG[i] = 1.0;
        }

        delcx.setup(natoms, coord, radii, coefS, coefV, coefM, coefG, vertices, tetra);

        delcx.regular3D(vertices, tetra, delaunay_eps);

        /*	==========================================================================================
            Generate alpha complex (with alpha=0.0)
            ========================================================================================== */

        double alpha = 0;
        alfcx.alfcx(alpha, vertices, tetra);

        /*	==========================================================================================
            Compute surface area and, optionally volume of the union of balls.
            If requested, compute also their derivatives
            ========================================================================================== */

        std::vector<Edge> edges;
        std::vector<Face> faces;
        alfcx.alphacxEdges(tetra, edges);
        alfcx.alphacxFaces(tetra, faces);

        double Surf, WSurf, Vol, WVol, Mean, WMean, Gauss, WGauss;

        int nfudge = 8;

        double *ballwvol, *dvol;
        ballwvol = new double[natoms+nfudge];
        dvol = new double[3*(natoms+nfudge)];
        memset(dvol, 0, 3*(natoms+nfudge)*sizeof(double));

        double *ballwsurf = new double[natoms+nfudge];
        double *dsurf = new double[3*(natoms+nfudge)];
        memset(dsurf, 0, 3*(natoms+nfudge)*sizeof(double));


        double *ballwmean, *dmean;
        ballwmean = new double[natoms+nfudge];
        dmean = new double[3*(natoms+nfudge)];
        memset(dmean, 0, 3*(natoms+nfudge)*sizeof(double));

        double *ballwgauss, *dgauss;
        ballwgauss = new double[natoms+nfudge];
        dgauss = new double[3*(natoms+nfudge)];
        memset(dgauss, 0, 3*(natoms+nfudge)*sizeof(double));

        volumes.ball_dvolumes(vertices, tetra, edges, faces, &WSurf, &WVol,
                              &WMean, &WGauss, &Surf, &Vol, &Mean, &Gauss, ballwsurf, ballwvol,
                              ballwmean, ballwgauss, dsurf, dvol, dmean, dgauss, 1);

        measures_out[0] = Vol;
        measures_out[1] = Surf;
        measures_out[2] = Mean;
        measures_out[3] = Gauss;

        for(int i = 0; i < natoms*3; i++){
            dvol_out[i] = dvol[i];
            dsurf_out[i] = dsurf[i];
            dmean_out[i] = dmean[i];
            dgauss_out[i] = dgauss[i];
        }

        delete [] coord; delete [] radii; delete [] coefS; delete [] coefV; delete [] coefM; delete [] coefG;
        delete [] ballwsurf; delete [] dsurf; delete [] ballwvol; delete [] dvol;
        delete [] ballwmean; delete [] dmean; delete [] ballwgauss; delete [] dgauss;
    });

    mod.method("get_geometric_measures_and_overlap_value_with_derivatives",  [](
            jlcxx::ArrayRef<double> measures_out,
            jlcxx::ArrayRef<double> dvol_out,
            jlcxx::ArrayRef<double> dsurf_out,
            jlcxx::ArrayRef<double> dmean_out,
            jlcxx::ArrayRef<double> dgauss_out,
            jlcxx::ArrayRef<double> dlol_out,
            const jlcxx::ArrayRef<double> in_coordinates, 
            const int molecule_size, 
            const jlcxx::ArrayRef<double> in_radii, 
            const double probe_radius,
            const double overlap_existence_penalty,
            const double overlap_penalty_slope,
            const double delaunay_eps){


        /*	==========================================================================================
            Compute Delaunay triangulation
            ========================================================================================== */

        clock_t start_s, stop_s;

        std::vector<Vertex> vertices;
        std::vector<Tetrahedron> tetra;

        int natoms = in_coordinates.size()/3;

        double *coord = new double[3*natoms];
        double *radii = new double[natoms];
        double *coefS = new double[natoms];
        double *coefV = new double[natoms];
        double *coefM = new double[natoms];
        double *coefG = new double[natoms];



        for(int i = 0; i < natoms; i++) {
            for(int j = 0; j < 3; j++){
                coord[3*i+j] = in_coordinates[3*i+j];
            }
            radii[i] = in_radii[i] + probe_radius;
            coefV[i] = 1.0;
            coefS[i] = 1.0;
            coefM[i] = 1.0;
            coefG[i] = 1.0;
        }

        delcx.setup(natoms, coord, radii, coefS, coefV, coefM, coefG, vertices, tetra);

        delcx.regular3D(vertices, tetra, delaunay_eps);

        /*	==========================================================================================
            Generate alpha complex (with alpha=0.0)
            ========================================================================================== */

        double alpha = 0;
        alfcx.alfcx(alpha, vertices, tetra);

        /*	==========================================================================================
            Compute surface area and, optionally volume of the union of balls.
            If requested, compute also their derivatives
            ========================================================================================== */

        std::vector<Edge> edges;
        std::vector<Face> faces;
        alfcx.alphacxEdges(tetra, edges);
        alfcx.alphacxFaces(tetra, faces);

        double Surf, WSurf, Vol, WVol, Mean, WMean, Gauss, WGauss;

        int nfudge = 8;

        double *ballwvol, *dvol;
        ballwvol = new double[natoms+nfudge];
        dvol = new double[3*(natoms+nfudge)];
        memset(dvol, 0, 3*(natoms+nfudge)*sizeof(double));

        double *ballwsurf = new double[natoms+nfudge];
        double *dsurf = new double[3*(natoms+nfudge)];
        memset(dsurf, 0, 3*(natoms+nfudge)*sizeof(double));


        double *ballwmean, *dmean;
        ballwmean = new double[natoms+nfudge];
        dmean = new double[3*(natoms+nfudge)];
        memset(dmean, 0, 3*(natoms+nfudge)*sizeof(double));

        double *ballwgauss, *dgauss;
        ballwgauss = new double[natoms+nfudge];
        dgauss = new double[3*(natoms+nfudge)];
        memset(dgauss, 0, 3*(natoms+nfudge)*sizeof(double));

        volumes.ball_dvolumes(vertices, tetra, edges, faces, &WSurf, &WVol,
                              &WMean, &WGauss, &Surf, &Vol, &Mean, &Gauss, ballwsurf, ballwvol,
                              ballwmean, ballwgauss, dsurf, dvol, dmean, dgauss, 1);

        double total_overlap_penalty = 0.0;
        int mol_size = 2;
        
        for (int i = 0; i < edges.size(); i++) {
            int ia = edges[i].Vertices[0] - 4;
            int ib = edges[i].Vertices[1] - 4;
            // Check if the edge is between two different molecules
            if (ia / molecule_size != ib / molecule_size) {
                double edge_ol = in_radii[ia] + in_radii[ib] - edges[i].Length;
                if (edge_ol > 0.0) {
                    total_overlap_penalty += overlap_existence_penalty + (overlap_penalty_slope * edge_ol);
                    // Compute the derivative of the overlap penalty
                    for (int j = 0; j < 3; j++) {
                        dlol_out[3*ia+j] = overlap_penalty_slope * (coord[3*ib+j] - coord[3*ia+j]);
                        dlol_out[3*ib+j] = overlap_penalty_slope * (coord[3*ia+j] - coord[3*ib+j]);
                    }
                }
            }
        }

        measures_out[0] = Vol;
        measures_out[1] = Surf;
        measures_out[2] = Mean;
        measures_out[3] = Gauss;
        measures_out[4] = total_overlap_penalty;

        for(int i = 0; i < natoms*3; i++){
            dvol_out[i] = dvol[i];
            dsurf_out[i] = dsurf[i];
            dmean_out[i] = dmean[i];
            dgauss_out[i] = dgauss[i];
        }

        delete [] coord; delete [] radii; delete [] coefS; delete [] coefV; delete [] coefM; delete [] coefG;
        delete [] ballwsurf; delete [] dsurf; delete [] ballwvol; delete [] dvol;
        delete [] ballwmean; delete [] dmean; delete [] ballwgauss; delete [] dgauss;
    });
}




