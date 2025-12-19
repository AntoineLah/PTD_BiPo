//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Tue Jul  8 17:19:42 2025 by ROOT version6.36.000)
//   from TTree particules/PTD data informations extracted ==> for second analysis
//   found on file: PTD_extracted_without_cut_for_elec_crossing_VD.root
//////////////////////////////////////////////////////////
#include <fstream>
#include <TF1.h>
#include <TGraph.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <limits>
#include <string>
#include <TH1F.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <climits>
//#include "/home/antoine/Documents/These/include_snome/sndisplay.cc"
//#include "/home/antoine/Documents/These/include_snome/mydict.cxx"
#include <fstream>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <filesystem>
using namespace std;
namespace fs = std::filesystem;
#include <TF1.h>

#include <fstream>
#include <TGraph.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <limits>
#include <string>
#include <TH1F.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <ctime>
#include <TPaveStats.h>
#include <TLegend.h>

#include "TFile.h"
#include "TTree.h"
#include "TLatex.h"
#include "TVector3.h"
#include "TSystem.h"
#include <sstream>
#include <iterator>
#include <regex>
#include <numeric>
#include <array>


inline std::vector<double> getRunNumbers(const std::string &folder) {
    std::regex pattern(R"(PTD_extracted_(\d+)\.root)");
    std::vector<double> runs;

    for (const auto &entry : fs::directory_iterator(folder)) {
        if (!entry.is_regular_file()) continue;

        std::string filename = entry.path().filename().string();
        std::smatch match;

        if (std::regex_match(filename, match, pattern)) {
            runs.push_back(std::stod(match[1].str())); // extract number
        }
    }

    return runs;
}


double calculateDistance(double x1, double y1, double z1, double x2, double y2, double z2) {
    return std::sqrt(std::pow(x1 - x2, 2) +
                     std::pow(y1 - y2, 2) +
                     std::pow(z1 - z2, 2));
}

struct RunInfo {
    bool found;         // true if trip found
    double tripTime;    // valid if found=true
    double run_start;   // always filled
    double duration;    // always filled
};

RunInfo getRunInfo(int run_number) {
    RunInfo info{false, 0.0, 0.0, 0.0};

    const char *run_list_filename;
    if (run_number >= 1546 && run_number < 2000)
        run_list_filename = "/sps/nemo/snemo/snemo_data/reco_data/UDD_betabeta_v1.list";
    else
        run_list_filename = "/sps/nemo/snemo/snemo_data/reco_data/UDD_betabeta_v2.list";

    std::ifstream run_list_file(run_list_filename);
    if (!run_list_file.is_open()) {
        std::cerr << "Error opening RunList file!" << std::endl;
        return info;
    }

    std::string header;
    std::getline(run_list_file, header); // skip header line

    int r;
    double run_start, duration, stop;
    std::string comment;

    while (run_list_file >> r >> run_start >> duration >> stop) {
        std::getline(run_list_file, comment); // read the rest of the line (comment text)

        if (r == run_number) {
            info.run_start = run_start;
            info.duration  = duration;

            if (stop != 0) {
                info.found    = true;
                info.tripTime = stop;
            }
            break;
        }
    }

    return info;
}



// help function to find the minimum distance between two tracks ! 
double find_min_distance_between_tracks(int elec_idx, int und_idx, 
	const vector<vector<double>>* x_start_per_elec_cluster,
	const vector<vector<double>>* y_start_per_elec_cluster,
	const vector<vector<double>>* z_start_per_elec_cluster,
	const vector<vector<double>>* x_end_per_elec_cluster,
	const vector<vector<double>>* y_end_per_elec_cluster,
	const vector<vector<double>>* z_end_per_elec_cluster,
	const vector<vector<double>>* x_start_per_UND_cluster,
	const vector<vector<double>>* y_start_per_UND_cluster,
	const vector<vector<double>>* z_start_per_UND_cluster,
	const vector<vector<double>>* x_end_per_UND_cluster,
	const vector<vector<double>>* y_end_per_UND_cluster,
	const vector<vector<double>>* z_end_per_UND_cluster,
	int elec_fit_idx,
	int und_fit_idx) 
{
		// Calculate distances between all possible combinations of endpoints
		double d1 = calculateDistance(
		x_start_per_elec_cluster->at(elec_idx).at(elec_fit_idx),
		y_start_per_elec_cluster->at(elec_idx).at(elec_fit_idx),
		z_start_per_elec_cluster->at(elec_idx).at(elec_fit_idx),
		x_start_per_UND_cluster->at(und_idx).at(und_fit_idx),
		y_start_per_UND_cluster->at(und_idx).at(und_fit_idx),
		z_start_per_UND_cluster->at(und_idx).at(und_fit_idx)
		);

		double d2 = calculateDistance(
		x_start_per_elec_cluster->at(elec_idx).at(elec_fit_idx),
		y_start_per_elec_cluster->at(elec_idx).at(elec_fit_idx),
		z_start_per_elec_cluster->at(elec_idx).at(elec_fit_idx),
		x_end_per_UND_cluster->at(und_idx).at(und_fit_idx),
		y_end_per_UND_cluster->at(und_idx).at(und_fit_idx),
		z_end_per_UND_cluster->at(und_idx).at(und_fit_idx)
		);

		double d3 = calculateDistance(
		x_end_per_elec_cluster->at(elec_idx).at(elec_fit_idx),
		y_end_per_elec_cluster->at(elec_idx).at(elec_fit_idx),
		z_end_per_elec_cluster->at(elec_idx).at(elec_fit_idx),
		x_start_per_UND_cluster->at(und_idx).at(und_fit_idx),
		y_start_per_UND_cluster->at(und_idx).at(und_fit_idx),
		z_start_per_UND_cluster->at(und_idx).at(und_fit_idx)
		);

		double d4 = calculateDistance(
		x_end_per_elec_cluster->at(elec_idx).at(elec_fit_idx),
		y_end_per_elec_cluster->at(elec_idx).at(elec_fit_idx),
		z_end_per_elec_cluster->at(elec_idx).at(elec_fit_idx),
		x_end_per_UND_cluster->at(und_idx).at(und_fit_idx),
		y_end_per_UND_cluster->at(und_idx).at(und_fit_idx),
		z_end_per_UND_cluster->at(und_idx).at(und_fit_idx)
		);

		// Return the minimum distance
		return std::min({d1, d2, d3, d4});
}

// Struct to represent a 3D vector 
struct Vec3 {
    double x, y, z;

    Vec3 operator-(const Vec3& v) const { return {x - v.x, y - v.y, z - v.z}; }
    Vec3 operator+(const Vec3& v) const { return {x + v.x, y + v.y, z + v.z}; }
    Vec3 operator*(double k) const { return {x * k, y * k, z * k}; }
};

// Scalar product
double dot(const Vec3& a, const Vec3& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

// Norme
double norm(const Vec3& v) {
    return sqrt(dot(v, v));
}
//Compute the closest points between two lines 
bool closestPointsBetweenLines(
    const Vec3& A, const Vec3& B,
    const Vec3& A2, const Vec3& B2,
    Vec3& C, Vec3& C2)
{
    Vec3 u = B - A;     // direction de la première droite
    Vec3 v = B2 - A2;   // direction de la deuxième droite
    Vec3 w0 = A - A2;   // vecteur entre les points de départ

    double a = dot(u, u);     // u·u
    double b = dot(u, v);     // u·v
    double c = dot(v, v);     // v·v
    double d = dot(u, w0);    // u·w0
    double e = dot(v, w0);    // v·w0
    double denom = a*c - b*b; // dénominateur

    // Si denom == 0 --> droites parallèles
    if (fabs(denom) < 1e-9) {
        return false; // parallèle ou quasi-parallèle
    }

    double t = (b*e - c*d) / denom;
    double s = (a*e - b*d) / denom;

    C  = A  + u * t;   // point sur D
    C2 = A2 + v * s;   // point sur D'

    return true;
}



void PTD_BiPo(int run){
        vector<double> OM2kill;    

        if(run >=1546 && run <= 1798 ){ // phase 0
            std::ifstream inputFile("/sps/nemo/scratch/granjon/full_gain_analysis/Bi/root/compute_a/phase_by_phase_calib/list_om_flag_phase_0.txt");  // open file
            if (!inputFile) {
                std::cerr << "Error opening file!" << std::endl;
            }
            
            double number;

            // Read doubles from file
            while (inputFile >> number) {
                OM2kill.push_back(number);
            }
            inputFile.close();
        }

        if(run >= 2000 && run <=2672){ // phase 1
            std::ifstream inputFile("/sps/nemo/scratch/granjon/full_gain_analysis/Bi/root/compute_a/phase_by_phase_calib/list_om_flag_phase_1.txt");  // open file
            if (!inputFile) {
                std::cerr << "Error opening file!" << std::endl;
            }
            
            double number;

            // Read doubles from file
            while (inputFile >> number) {
                OM2kill.push_back(number);
            }
            inputFile.close();
        }
        if(run >= 2673 && run <= 2869){ // phase 2a
            std::ifstream inputFile("/sps/nemo/scratch/granjon/full_gain_analysis/Bi/root/compute_a/phase_by_phase_calib/list_om_phase_2a.txt");  // open file
            if (!inputFile) {
                std::cerr << "Error opening file!" << std::endl;
            }
            std::string line;
            std::getline(inputFile, line);
            
            double number;

            // Read doubles from file
            while (inputFile >> number) {
                OM2kill.push_back(number);
            }
            inputFile.close();
        }

        if(run >= 2870 && run <= 3467){ // phase 2b 
            std::ifstream inputFile("/sps/nemo/scratch/granjon/full_gain_analysis/Bi/root/compute_a/phase_by_phase_calib/list_om_phase_2b.txt");  // open file
            if (!inputFile) {
                std::cerr << "Error opening file!" << std::endl;
            }
            std::string line;
            std::getline(inputFile, line);
            
            double number;

            // Read doubles from file
            while (inputFile >> number) {
                OM2kill.push_back(number);
            }
            inputFile.close();

        }
        if(run >= 3470 && run < 5000){ // phase 3
            std::ifstream inputFile("/sps/nemo/scratch/granjon/full_gain_analysis/Bi/root/compute_a/phase_by_phase_calib/list_om_phase_3.txt");  // open file
            if (!inputFile) {
                std::cerr << "Error opening file!" << std::endl;
            }
            
            double number;

            // Read doubles from file
            while (inputFile >> number) {
                OM2kill.push_back(number);
            }
            inputFile.close();
        }
        std::string folder ="/sps/nemo/scratch/lahaie/FalaiseSkeletonModules/build/out/renamed_files/";
    
        std::cout<<"run ="<<run<<std::endl;
        std::string file_path = folder + "PTD_extracted_" + std::to_string((int)run) + ".root";
        std::cout << "Processing file: " << file_path << std::endl;
        
        RunInfo info = getRunInfo(run);
        
        //sndisplay::demonstrator *sndemonstrator = new sndisplay::demonstrator ("demonstrator_test");
        //sndemonstrator->setrange(0,1);
        // Open the ROOT file
        TFile *f = TFile::Open(file_path.c_str(), "READ");
        if(!f || f->IsZombie()){
            std::cerr << "Error: cannot open input file " << file_path << std::endl;
            if(f){ f->Close(); delete f; }
            
        }

        std::string outFileName = "/sps/nemo/scratch/lahaie/FalaiseSkeletonModules/build/out/BiBo/all_phases/" + std::to_string((int)run) + "_out.root";
        //std::string outFileName = "/sps/nemo/scratch/lahaie/FalaiseSkeletonModules/build/out/BiBo/injection/" + std::to_string((int)run) + "_out.root";
        // ensure output directory exists
        try {
            fs::path p(outFileName);
            if(!p.has_parent_path()){
                std::cerr << "Warning: output path has no parent: " << outFileName << std::endl;
            } else {
                fs::create_directories(p.parent_path());
            }
        } catch(const std::exception &e){
            std::cerr << "Error creating output directory: " << e.what() << std::endl;
            f->Close(); delete f;
            
        }

        TFile *outFile = TFile::Open(outFileName.c_str(), "RECREATE");
        if(!outFile || outFile->IsZombie()){
            std::cerr << "Error: cannot create output file " << outFileName << std::endl;
            if(outFile){ outFile->Close(); delete outFile; }
            f->Close(); delete f;
            
        }

        // Get the TTree
        TTree *particules = dynamic_cast<TTree*>(f->Get("particules"));
        if (!particules) {
           std::cerr << "Error: Cannot find the TTree 'particules' in the file " << file_path << " -- skipping file." << std::endl;
           outFile->Close(); delete outFile;
           f->Close();
           
        }

        

   
        //Declaration of leaves types
        int event_number ;
        int nb_isolated_calo ;
        vector<double>  *E_isolated_calo = nullptr;
        vector<double>  *time_isolated_calo = nullptr;
        vector<double>  *sigma_E_isolated_calo = nullptr;
        vector<double>  *sigma_time_isolated_calo = nullptr;
        vector<char>    *isolated_calo_type = nullptr;
        vector<double>  *isolated_calo_timestamp = nullptr;
        vector<double>  *isolated_calo_charge = nullptr;
        vector<double>  *isolated_calo_amplitude = nullptr;
        vector<double>  *isolated_calo_num = nullptr;
        vector<double>  *sigma_isolated_calo_timestamp = nullptr;
        vector<double>  *sigma_isolated_calo_charge = nullptr;
        vector<double>  *sigma_isolated_calo_amplitude = nullptr;
        bool          evt_isolated_calo;
        vector<bool>    *isolated_calo_low_threshold_only = nullptr;
        vector<bool>    *isolated_calo_high_threshold = nullptr;
        vector<double> *x_calo_gamma = nullptr; 
        vector<double> *y_calo_gamma = nullptr;
        vector<double> *z_calo_gamma = nullptr; 
        bool          tracks_with_associated_calo ;
        int           nb_of_elec_candidates;
        vector<vector<double> > *cell_num_per_elec_cluster = nullptr;
        vector<vector<double> > *OM_num_per_elec_cluster = nullptr;
        vector<vector<double> > *E_OM_per_elec_cluster = nullptr;
        vector<double>           *nb_of_OM_per_elec_cluster = nullptr;
        vector<vector<double> > *anode_time_per_elec_cluster = nullptr;
        vector<vector<double> > *top_cathode_per_elec_cluster = nullptr;
        vector<vector<double> > *bottom_cathode_per_elec_cluster = nullptr;
        vector<vector<double> > *OM_timestamp_per_elec_cluster = nullptr;
        vector<vector<double> > *OM_charge_per_elec_cluster = nullptr;
        vector<vector<double> > *OM_amplitude_per_elec_cluster = nullptr;
        vector<vector<double> > *OM_LT_only_per_elec_cluster = nullptr;
        vector<vector<double> > *OM_HT_per_elec_cluster= nullptr;
        vector<vector<double> > *z_of_cells_per_elec_cluster= nullptr;
        vector<vector<double> > *sigma_z_of_cells_per_elec_cluster= nullptr;
        vector<vector<double> > *r_of_cells_per_elec_cluster= nullptr;
        vector<vector<double> > *sigma_r_of_cells_per_elec_cluster= nullptr;
        vector<bool>    *elec_cluster_is_delayed = nullptr;
        vector<bool>    *elec_clsuter_is_prompt = nullptr;
        vector<vector<double> > *x_start_per_elec_cluster = nullptr;
        vector<vector<double> > *y_start_per_elec_cluster = nullptr;
        vector<vector<double> > *z_start_per_elec_cluster = nullptr;
        vector<vector<double> > *x_end_per_elec_cluster = nullptr;
        vector<vector<double> > *y_end_per_elec_cluster = nullptr;
        vector<vector<double> > *z_end_per_elec_cluster = nullptr;
        vector<int>     *ID_clsuter_per_elec_cluster = nullptr;
        vector<int>     *nb_fit_solution_per_elec_cluster = nullptr;
        vector<double>  *nb_cell_per_elec_cluster = nullptr;
        vector<vector<double> > *x_is_on_reference_source_plane = nullptr;
        vector<vector<double> > *y_is_on_reference_source_plane = nullptr;
        vector<vector<double> > *z_is_on_reference_source_plane = nullptr;
        vector<vector<double> > *x_is_on_source_foil = nullptr;
        vector<vector<double> > *y_is_on_source_foil = nullptr;
        vector<vector<double> > *z_is_on_source_foil = nullptr;
        vector<vector<double> > *x_is_on_main_calorimeter = nullptr;
        vector<vector<double> > *y_is_on_main_calorimeter = nullptr;
        vector<vector<double> > *z_is_on_main_calorimeter = nullptr;
        vector<vector<double> > *x_is_on_x_calorimeter = nullptr;
        vector<vector<double> > *y_is_on_x_calorimeter = nullptr;
        vector<vector<double> > *z_is_on_x_calorimeter = nullptr;
        vector<vector<double> > *x_is_on_gamma_veto = nullptr;
        vector<vector<double> > *y_is_on_gamma_veto = nullptr;
        vector<vector<double> > *z_is_on_gamma_veto = nullptr;
        vector<vector<double> > *x_is_on_wire = nullptr;
        vector<vector<double> > *y_is_on_wire = nullptr;
        vector<vector<double> > *z_is_on_wire = nullptr;
        vector<vector<double> > *x_is_in_gas = nullptr;
        vector<vector<double> > *y_is_in_gas = nullptr;
        vector<vector<double> > *z_is_in_gas = nullptr;
        vector<vector<double> > *x_is_on_source_gap = nullptr;
        vector<vector<double> > *y_is_on_source_gap = nullptr;
        vector<vector<double> > *z_is_on_source_gap = nullptr;
        vector<vector<double> > *xy_distance_from_reference_source_plane = nullptr;
        vector<vector<double> > *xy_distance_from_source_foil = nullptr;
        vector<vector<double> > *xy_distance_from_main_calorimeter = nullptr;
        vector<vector<double> > *xy_distance_from_x_calorimeter = nullptr;
        vector<vector<double> > *xy_distance_from_gamma_veto = nullptr;
        vector<vector<double> > *xy_distance_from_wire = nullptr;
        vector<vector<double> > *xy_distance_from_gas = nullptr;
        vector<vector<double> > *xy_distance_from_source_gap = nullptr;
        vector<vector<double> > *xyz_distance_from_reference_source_plane = nullptr;
        vector<vector<double> > *xyz_distance_from_source_foil = nullptr;
        vector<vector<double> > *xyz_distance_from_main_calorimeter = nullptr;
        vector<vector<double> > *xyz_distance_from_x_calorimeter = nullptr;
        vector<vector<double> > *xyz_distance_from_gamma_veto = nullptr;
        vector<vector<double> > *xyz_distance_from_wire = nullptr;
        vector<vector<double> > *xyz_distance_from_gas = nullptr;
        vector<vector<double> > *xyz_distance_from_source_gap = nullptr;
        vector<vector<double> > *vertex_is_on_reference_source_plane_per_elec_cluster = nullptr;
        vector<vector<double> > *vertex_is_on_source_foil_per_elec_cluster = nullptr;
        vector<vector<double> > *vertex_is_on_main_calorimeter_per_elec_cluster = nullptr;
        vector<vector<double> > *vertex_is_on_x_calorimeter_per_elec_cluster = nullptr;
        vector<vector<double> > *vertex_is_on_gamma_veto_per_elec_cluster = nullptr;
        vector<vector<double> > *vertex_is_on_wire_per_elec_cluster = nullptr;
        vector<vector<double> > *vertex_is_in_gas_per_elec_cluster = nullptr;
        vector<vector<double> > *vertex_is_on_source_gap_per_elec_cluster = nullptr;
        vector<vector<double> > *distance_elec_UND_per_elec_crossing = nullptr;
        int           nb_elec_crossing;
        vector<vector<double> > *delta_t_cells_of_UND_per_elec_crossing = nullptr;
        vector<vector<double> > *delta_t_cells_of_elec_per_elec_crossing = nullptr;
        vector<bool>    *elec_used = nullptr;
        vector<bool>    *und_used = nullptr;
        vector<double>  *nb_of_cell_per_elec_crossing = nullptr;
        vector<int>     *ID_cluster_per_elec_crossing = nullptr;
        vector<int>     *ID_cluster_UND_per_elec_crossing = nullptr;
        bool          tracks_without_associated_calo;

        vector<vector<double>> *und_elec_crossing_xs = nullptr;
	    vector<vector<double>> *und_elec_crossing_ys = nullptr;
	    vector<vector<double>> *und_elec_crossing_zs = nullptr;

	    vector<vector<double>> *und_elec_crossing_xe = nullptr ;
	    vector<vector<double>> *und_elec_crossing_ye = nullptr ;
	    vector<vector<double>> *und_elec_crossing_ze = nullptr ;

        int           nb_of_UND_candidates;
        vector<vector<double> > *cell_num_per_UND_cluster = nullptr;
        vector<vector<double> > *anode_time_per_UND_cluster = nullptr;
        vector<vector<double> > *top_cathodes_per_UND_cluster = nullptr;
        vector<vector<double> > *bottom_cathodes_per_UND_cluster = nullptr;
        vector<vector<double> > *z_of_cell_per_UND_cluster = nullptr;
        vector<vector<double> > *sigma_z_of_cells_per_UND_cluster = nullptr;
        vector<vector<double> > *r_of_cells_per_UND_cluster = nullptr ;
        vector<vector<double> > *sigma_r_of_cells_per_UND_cluster = nullptr;
        vector<bool>    *UND_cluster_is_delayed = nullptr;
        vector<bool>    *UND_cluster_is_prompt = nullptr;
        vector<int>     *ID_clsuter_UND = nullptr;
        vector<vector<double> > *x_start_per_UND_cluster = nullptr;
        vector<vector<double> > *y_start_per_UND_cluster = nullptr;
        vector<vector<double> > *z_start_per_UND_cluster = nullptr;
        vector<vector<double> > *x_end_per_UND_cluster = nullptr;
        vector<vector<double> > *y_end_per_UND_cluster = nullptr;
        vector<vector<double> > *z_end_per_UND_cluster = nullptr;
        vector<int>     *nb_cell_per_UND_cluster = nullptr;
        vector<double>  *min_time_per_UND_cluster = nullptr;
        vector<int>     *nb_fit_solution_per_UND_cluster = nullptr;
        vector<vector<double> > *x_is_on_reference_source_plane_per_UND_cluster = nullptr;
        vector<vector<double> > *y_is_on_reference_source_plane_per_UND_cluster = nullptr;
        vector<vector<double> > *z_is_on_reference_source_plane_per_UND_cluster = nullptr;
        vector<vector<double> > *x_is_on_source_foil_per_UND_cluster = nullptr;
        vector<vector<double> > *y_is_on_source_foil_per_UND_cluster = nullptr;
        vector<vector<double> > *z_is_on_source_foil_per_UND_cluster = nullptr;
        vector<vector<double> > *x_is_on_wire_per_UND_cluster = nullptr;
        vector<vector<double> > *y_is_on_wire_per_UND_cluster = nullptr;
        vector<vector<double> > *z_is_on_wire_per_UND_cluster = nullptr;
        vector<vector<double> > *x_is_in_gas_per_UND_cluster = nullptr;
        vector<vector<double> > *y_is_in_gas_per_UND_cluster = nullptr;
        vector<vector<double> > *z_is_in_gas_per_UND_cluster = nullptr;
        vector<vector<double> > *x_is_on_source_gap_per_UND_cluster = nullptr;
        vector<vector<double> > *y_is_on_source_gap_per_UND_cluster = nullptr;
        vector<vector<double> > *z_is_on_source_gap_per_UND_cluster = nullptr;
        vector<vector<double> > *xy_distance_from_reference_source_plane_per_UND_cluster = nullptr;
        vector<vector<double> > *xy_distance_from_source_foil_per_UND_cluster = nullptr;
        vector<vector<double> > *xy_distance_from_wire_per_UND_cluster = nullptr;
        vector<vector<double> > *xy_distance_from_gas_per_UND_cluster = nullptr;
        vector<vector<double> > *xy_distance_from_source_gap_per_UND_cluster = nullptr;
        vector<vector<double> > *xyz_distance_from_reference_source_plane_per_UND_cluster = nullptr;
        vector<vector<double> > *xyz_distance_from_source_foil_per_UND_cluster = nullptr;
        vector<vector<double> > *xyz_distance_from_wire_per_UND_cluster = nullptr;
        vector<vector<double> > *xyz_distance_from_gas_per_UND_cluster = nullptr;
        vector<vector<double> > *xyz_distance_from_source_gap_per_UND_cluster = nullptr;
        bool          UND_is_on_reference_source_plane;
        bool          UND_is_on_source_foil;
        bool          UND_is_on_wire;
        bool          UND_is_in_gas;
        bool          UND_is_on_source_gap;
        vector<vector<double> > *vertex_is_on_reference_source_plane_per_UND_cluster = nullptr;
        vector<vector<double> > *vertex_is_on_source_foil_per_UND_cluster = nullptr;
        vector<vector<double> > *vertex_is_on_wire_per_UND_cluster = nullptr;
        vector<vector<double> > *vertex_is_in_gas_per_UND_cluster = nullptr;
        vector<vector<double> > *vertex_is_on_source_gap_per_UND_cluster = nullptr;
        vector<int> *elec_is_associated_with_track = nullptr;
 


        double couleur;

        // Set branch addresses.
        particules->SetBranchAddress("event_number",&event_number);
        particules->SetBranchAddress("nb_isolated_calo",&nb_isolated_calo);
        particules->SetBranchAddress("E_isolated_calo",&E_isolated_calo);
        particules->SetBranchAddress("time_isolated_calo",&time_isolated_calo);
        particules->SetBranchAddress("sigma_E_isolated_calo",&sigma_E_isolated_calo);
        particules->SetBranchAddress("sigma_time_isolated_calo",&sigma_time_isolated_calo);
        particules->SetBranchAddress("isolated_calo_type",&isolated_calo_type);
        particules->SetBranchAddress("isolated_calo_timestamp",&isolated_calo_timestamp);
        particules->SetBranchAddress("isolated_calo_charge",&isolated_calo_charge);
        particules->SetBranchAddress("isolated_calo_amplitude",&isolated_calo_amplitude);
        particules->SetBranchAddress("isolated_calo_num",&isolated_calo_num);
        particules->SetBranchAddress("sigma_isolated_calo_timestamp",&sigma_isolated_calo_timestamp);
        particules->SetBranchAddress("sigma_isolated_calo_charge",&sigma_isolated_calo_charge);
        particules->SetBranchAddress("sigma_isolated_calo_amplitude",&sigma_isolated_calo_amplitude);
        particules->SetBranchAddress("evt_isolated_calo",&evt_isolated_calo);
        particules->SetBranchAddress("isolated_calo_low_threshold_only",&isolated_calo_low_threshold_only);
        particules->SetBranchAddress("isolated_calo_high_threshold",&isolated_calo_high_threshold);
        particules->SetBranchAddress("x_calo_gamma",&x_calo_gamma);
        particules->SetBranchAddress("y_calo_gamma",&y_calo_gamma);
        particules->SetBranchAddress("z_calo_gamma",&z_calo_gamma);

        particules->SetBranchAddress("tracks_with_associated_calo",&tracks_with_associated_calo);
        particules->SetBranchAddress("nb_of_elec_candidates",&nb_of_elec_candidates);
        particules->SetBranchAddress("cell_num_per_elec_cluster",&cell_num_per_elec_cluster);
        particules->SetBranchAddress("OM_num_per_elec_cluster",&OM_num_per_elec_cluster);
        particules->SetBranchAddress("E_OM_per_elec_cluster",&E_OM_per_elec_cluster);
        particules->SetBranchAddress("nb_of_OM_per_elec_cluster",&nb_of_OM_per_elec_cluster);
        particules->SetBranchAddress("anode_time_per_elec_cluster",&anode_time_per_elec_cluster);
        particules->SetBranchAddress("top_cathode_per_elec_cluster",&top_cathode_per_elec_cluster);
        particules->SetBranchAddress("bottom_cathode_per_elec_cluster",&bottom_cathode_per_elec_cluster);
        particules->SetBranchAddress("OM_timestamp_per_elec_cluster",&OM_timestamp_per_elec_cluster);
        particules->SetBranchAddress("OM_charge_per_elec_cluster",&OM_charge_per_elec_cluster);
        particules->SetBranchAddress("OM_amplitude_per_elec_cluster",&OM_amplitude_per_elec_cluster);
        particules->SetBranchAddress("OM_LT_only_per_elec_cluster",&OM_LT_only_per_elec_cluster);
        particules->SetBranchAddress("OM_HT_per_elec_cluster",&OM_HT_per_elec_cluster);
        particules->SetBranchAddress("z_of_cells_per_elec_cluster",&z_of_cells_per_elec_cluster);
        particules->SetBranchAddress("sigma_z_of_cells_per_elec_cluster",&sigma_z_of_cells_per_elec_cluster);
        particules->SetBranchAddress("r_of_cells_per_elec_cluster",&r_of_cells_per_elec_cluster);
        particules->SetBranchAddress("sigma_r_of_cells_per_elec_cluster",&sigma_r_of_cells_per_elec_cluster);
        particules->SetBranchAddress("elec_cluster_is_delayed",&elec_cluster_is_delayed);
        particules->SetBranchAddress("elec_clsuter_is_prompt",&elec_clsuter_is_prompt);
        particules->SetBranchAddress("x_start_per_elec_cluster",&x_start_per_elec_cluster);
        particules->SetBranchAddress("y_start_per_elec_cluster",&y_start_per_elec_cluster);
        particules->SetBranchAddress("z_start_per_elec_cluster",&z_start_per_elec_cluster);
        particules->SetBranchAddress("x_end_per_elec_cluster",&x_end_per_elec_cluster);
        particules->SetBranchAddress("y_end_per_elec_cluster",&y_end_per_elec_cluster);
        particules->SetBranchAddress("z_end_per_elec_cluster",&z_end_per_elec_cluster);
        particules->SetBranchAddress("ID_clsuter_per_elec_cluster",&ID_clsuter_per_elec_cluster);

        particules->SetBranchAddress("nb_fit_solution_per_elec_cluster",&nb_fit_solution_per_elec_cluster);
        particules->SetBranchAddress("nb_cell_per_elec_cluster",&nb_cell_per_elec_cluster);
        particules->SetBranchAddress("x_is_on_reference_source_plane",&x_is_on_reference_source_plane);
        particules->SetBranchAddress("y_is_on_reference_source_plane",&y_is_on_reference_source_plane);
        particules->SetBranchAddress("z_is_on_reference_source_plane",&z_is_on_reference_source_plane);
        particules->SetBranchAddress("x_is_on_source_foil",&x_is_on_source_foil);
        particules->SetBranchAddress("y_is_on_source_foil",&y_is_on_source_foil);
        particules->SetBranchAddress("z_is_on_source_foil",&z_is_on_source_foil);
        particules->SetBranchAddress("x_is_on_main_calorimeter",&x_is_on_main_calorimeter);
        particules->SetBranchAddress("y_is_on_main_calorimeter",&y_is_on_main_calorimeter);
        particules->SetBranchAddress("z_is_on_main_calorimeter",&z_is_on_main_calorimeter);
        particules->SetBranchAddress("x_is_on_x_calorimeter",&x_is_on_x_calorimeter);
        particules->SetBranchAddress("y_is_on_x_calorimeter",&y_is_on_x_calorimeter);
        particules->SetBranchAddress("z_is_on_x_calorimeter",&z_is_on_x_calorimeter);
        particules->SetBranchAddress("x_is_on_gamma_veto",&x_is_on_gamma_veto);
        particules->SetBranchAddress("y_is_on_gamma_veto",&y_is_on_gamma_veto);
        particules->SetBranchAddress("z_is_on_gamma_veto",&z_is_on_gamma_veto);
        particules->SetBranchAddress("x_is_on_wire",&x_is_on_wire);
        particules->SetBranchAddress("y_is_on_wire",&y_is_on_wire);
        particules->SetBranchAddress("z_is_on_wire",&z_is_on_wire);
        particules->SetBranchAddress("x_is_in_gas",&x_is_in_gas);
        particules->SetBranchAddress("y_is_in_gas",&y_is_in_gas);
        particules->SetBranchAddress("z_is_in_gas",&z_is_in_gas);
        particules->SetBranchAddress("x_is_on_source_gap",&x_is_on_source_gap);
        particules->SetBranchAddress("y_is_on_source_gap",&y_is_on_source_gap);
        particules->SetBranchAddress("z_is_on_source_gap",&z_is_on_source_gap);
        particules->SetBranchAddress("xy_distance_from_reference_source_plane",&xy_distance_from_reference_source_plane);
        particules->SetBranchAddress("xy_distance_from_source_foil",&xy_distance_from_source_foil);
        particules->SetBranchAddress("xy_distance_from_main_calorimeter",&xy_distance_from_main_calorimeter);
        particules->SetBranchAddress("xy_distance_from_x_calorimeter",&xy_distance_from_x_calorimeter);
        particules->SetBranchAddress("xy_distance_from_gamma_veto",&xy_distance_from_gamma_veto);
        particules->SetBranchAddress("xy_distance_from_wire",&xy_distance_from_wire);
        particules->SetBranchAddress("xy_distance_from_gas",&xy_distance_from_gas);
        particules->SetBranchAddress("xy_distance_from_source_gap",&xy_distance_from_source_gap);
        particules->SetBranchAddress("xyz_distance_from_reference_source_plane",&xyz_distance_from_reference_source_plane);
        particules->SetBranchAddress("xyz_distance_from_source_foil",&xyz_distance_from_source_foil);
        particules->SetBranchAddress("xyz_distance_from_main_calorimeter",&xyz_distance_from_main_calorimeter);
        particules->SetBranchAddress("xyz_distance_from_x_calorimeter",&xyz_distance_from_x_calorimeter);
        particules->SetBranchAddress("xyz_distance_from_gamma_veto",&xyz_distance_from_gamma_veto);
        particules->SetBranchAddress("xyz_distance_from_wire",&xyz_distance_from_wire);
        particules->SetBranchAddress("xyz_distance_from_gas",&xyz_distance_from_gas);
        particules->SetBranchAddress("xyz_distance_from_source_gap",&xyz_distance_from_source_gap);
        particules->SetBranchAddress("vertex_is_on_reference_source_plane_per_elec_cluster",&vertex_is_on_reference_source_plane_per_elec_cluster);
        particules->SetBranchAddress("vertex_is_on_source_foil_per_elec_cluster",&vertex_is_on_source_foil_per_elec_cluster);
        particules->SetBranchAddress("vertex_is_on_main_calorimeter_per_elec_cluster",&vertex_is_on_main_calorimeter_per_elec_cluster);
        particules->SetBranchAddress("vertex_is_on_x_calorimeter_per_elec_cluster",&vertex_is_on_x_calorimeter_per_elec_cluster);
        particules->SetBranchAddress("vertex_is_on_gamma_veto_per_elec_cluster",&vertex_is_on_gamma_veto_per_elec_cluster);
        particules->SetBranchAddress("vertex_is_on_wire_per_elec_cluster",&vertex_is_on_wire_per_elec_cluster);
        particules->SetBranchAddress("vertex_is_in_gas_per_elec_cluster",&vertex_is_in_gas_per_elec_cluster);
        particules->SetBranchAddress("vertex_is_on_source_gap_per_elec_cluster",&vertex_is_on_source_gap_per_elec_cluster);

        particules->SetBranchAddress("distance_elec_UND_per_elec_crossing",&distance_elec_UND_per_elec_crossing);
        particules->SetBranchAddress("nb_elec_crossing",&nb_elec_crossing);
        particules->SetBranchAddress("delta_t_cells_of_UND_per_elec_crossing",&delta_t_cells_of_UND_per_elec_crossing);
        particules->SetBranchAddress("delta_t_cells_of_elec_per_elec_crossing",&delta_t_cells_of_elec_per_elec_crossing);

        particules->SetBranchAddress("elec_used",&elec_used);

        particules->SetBranchAddress("und_used",&und_used);
        particules->SetBranchAddress("nb_of_cell_per_elec_crossing",&nb_of_cell_per_elec_crossing);
        particules->SetBranchAddress("ID_cluster_per_elec_crossing",&ID_cluster_per_elec_crossing);
        particules->SetBranchAddress("ID_cluster_UND_per_elec_crossing",&ID_cluster_UND_per_elec_crossing);
        particules->SetBranchAddress("elec_is_associated_with_track", &elec_is_associated_with_track);

        particules->SetBranchAddress("und_elec_crossing_xs", &und_elec_crossing_xs);
        particules->SetBranchAddress("und_elec_crossing_ys", &und_elec_crossing_ys);
        particules->SetBranchAddress("und_elec_crossing_zs", &und_elec_crossing_zs);
        particules->SetBranchAddress("und_elec_crossing_xe", &und_elec_crossing_xe);
        particules->SetBranchAddress("und_elec_crossing_ye", &und_elec_crossing_ye);
        particules->SetBranchAddress("und_elec_crossing_ze", &und_elec_crossing_ze);

        particules->SetBranchAddress("tracks_without_associated_calo",&tracks_without_associated_calo);
        particules->SetBranchAddress("nb_of_UND_candidates",&nb_of_UND_candidates);
        particules->SetBranchAddress("cell_num_per_UND_cluster",&cell_num_per_UND_cluster);
        particules->SetBranchAddress("anode_time_per_UND_cluster",&anode_time_per_UND_cluster);
        particules->SetBranchAddress("top_cathodes_per_UND_cluster",&top_cathodes_per_UND_cluster);
        particules->SetBranchAddress("bottom_cathodes_per_UND_cluster",&bottom_cathodes_per_UND_cluster);
        particules->SetBranchAddress("z_of_cell_per_UND_cluster",&z_of_cell_per_UND_cluster);
        particules->SetBranchAddress("sigma_z_of_cells_per_UND_cluster",&sigma_z_of_cells_per_UND_cluster);
        particules->SetBranchAddress("r_of_cells_per_UND_cluster",&r_of_cells_per_UND_cluster);
        particules->SetBranchAddress("sigma_r_of_cells_per_UND_cluster",&sigma_r_of_cells_per_UND_cluster);
        particules->SetBranchAddress("UND_cluster_is_delayed",&UND_cluster_is_delayed);
        particules->SetBranchAddress("UND_cluster_is_prompt",&UND_cluster_is_prompt);
        particules->SetBranchAddress("ID_clsuter_UND",&ID_clsuter_UND);
        particules->SetBranchAddress("x_start_per_UND_cluster",&x_start_per_UND_cluster);
        particules->SetBranchAddress("y_start_per_UND_cluster",&y_start_per_UND_cluster);
        particules->SetBranchAddress("z_start_per_UND_cluster",&z_start_per_UND_cluster);
        particules->SetBranchAddress("x_end_per_UND_cluster",&x_end_per_UND_cluster);
        particules->SetBranchAddress("y_end_per_UND_cluster",&y_end_per_UND_cluster);
        particules->SetBranchAddress("z_end_per_UND_cluster",&z_end_per_UND_cluster);
        particules->SetBranchAddress("nb_cell_per_UND_cluster",&nb_cell_per_UND_cluster);
        particules->SetBranchAddress("min_time_per_UND_cluster",&min_time_per_UND_cluster);
        particules->SetBranchAddress("nb_fit_solution_per_UND_cluster",&nb_fit_solution_per_UND_cluster);
        particules->SetBranchAddress("x_is_on_reference_source_plane_per_UND_cluster",&x_is_on_reference_source_plane_per_UND_cluster);
        particules->SetBranchAddress("y_is_on_reference_source_plane_per_UND_cluster",&y_is_on_reference_source_plane_per_UND_cluster);
        particules->SetBranchAddress("z_is_on_reference_source_plane_per_UND_cluster",&z_is_on_reference_source_plane_per_UND_cluster);

        particules->SetBranchAddress("x_is_on_source_foil_per_UND_cluster",&x_is_on_source_foil_per_UND_cluster);
        particules->SetBranchAddress("y_is_on_source_foil_per_UND_cluster",&y_is_on_source_foil_per_UND_cluster);
        particules->SetBranchAddress("z_is_on_source_foil_per_UND_cluster",&z_is_on_source_foil_per_UND_cluster);
        particules->SetBranchAddress("x_is_on_wire_per_UND_cluster",&x_is_on_wire_per_UND_cluster);
        particules->SetBranchAddress("y_is_on_wire_per_UND_cluster",&y_is_on_wire_per_UND_cluster);
        particules->SetBranchAddress("z_is_on_wire_per_UND_cluster",&z_is_on_wire_per_UND_cluster);
        particules->SetBranchAddress("x_is_in_gas_per_UND_cluster",&x_is_in_gas_per_UND_cluster);
        particules->SetBranchAddress("y_is_in_gas_per_UND_cluster",&y_is_in_gas_per_UND_cluster);
        particules->SetBranchAddress("z_is_in_gas_per_UND_cluster",&z_is_in_gas_per_UND_cluster);
        particules->SetBranchAddress("x_is_on_source_gap_per_UND_cluster",&x_is_on_source_gap_per_UND_cluster);
        particules->SetBranchAddress("y_is_on_source_gap_per_UND_cluster",&y_is_on_source_gap_per_UND_cluster);
        particules->SetBranchAddress("z_is_on_source_gap_per_UND_cluster",&z_is_on_source_gap_per_UND_cluster);
        particules->SetBranchAddress("xy_distance_from_reference_source_plane_per_UND_cluster",&xy_distance_from_reference_source_plane_per_UND_cluster);
        particules->SetBranchAddress("xy_distance_from_source_foil_per_UND_cluster",&xy_distance_from_source_foil_per_UND_cluster);
        particules->SetBranchAddress("xy_distance_from_wire_per_UND_cluster",&xy_distance_from_wire_per_UND_cluster);
        particules->SetBranchAddress("xy_distance_from_gas_per_UND_cluster",&xy_distance_from_gas_per_UND_cluster);
        particules->SetBranchAddress("xy_distance_from_source_gap_per_UND_cluster",&xy_distance_from_source_gap_per_UND_cluster);
        particules->SetBranchAddress("xyz_distance_from_reference_source_plane_per_UND_cluster",&xyz_distance_from_reference_source_plane_per_UND_cluster);
        particules->SetBranchAddress("xyz_distance_from_source_foil_per_UND_cluster",&xyz_distance_from_source_foil_per_UND_cluster);
        particules->SetBranchAddress("xyz_distance_from_wire_per_UND_cluster",&xyz_distance_from_wire_per_UND_cluster);
        particules->SetBranchAddress("xyz_distance_from_gas_per_UND_cluster",&xyz_distance_from_gas_per_UND_cluster);
        particules->SetBranchAddress("xyz_distance_from_source_gap_per_UND_cluster",&xyz_distance_from_source_gap_per_UND_cluster);
        particules->SetBranchAddress("UND_is_on_reference_source_plane",&UND_is_on_reference_source_plane);
        particules->SetBranchAddress("UND_is_on_source_foil",&UND_is_on_source_foil);
        particules->SetBranchAddress("UND_is_on_wire",&UND_is_on_wire);
        particules->SetBranchAddress("UND_is_in_gas",&UND_is_in_gas);
        particules->SetBranchAddress("UND_is_on_source_gap",&UND_is_on_source_gap);
        particules->SetBranchAddress("vertex_is_on_reference_source_plane_per_UND_cluster",&vertex_is_on_reference_source_plane_per_UND_cluster);
        particules->SetBranchAddress("vertex_is_on_source_foil_per_UND_cluster",&vertex_is_on_source_foil_per_UND_cluster);
        particules->SetBranchAddress("vertex_is_on_wire_per_UND_cluster",&vertex_is_on_wire_per_UND_cluster);
        particules->SetBranchAddress("vertex_is_in_gas_per_UND_cluster",&vertex_is_in_gas_per_UND_cluster);
        particules->SetBranchAddress("vertex_is_on_source_gap_per_UND_cluster",&vertex_is_on_source_gap_per_UND_cluster);

        // create output ROOT file and TTree to store electron-UND pair information
    outFile->cd();
    TTree *pairTree = new TTree("pairs","Electron-UND fit pairs");

    // tree variables
    int t_event = 0;
    int t_elecIdx = -1, t_elecFit = -1, t_undIdx = -1, t_undFit = -1;
    int t_elecClusterID = -1, t_undClusterID = -1;
    double t_minDist = -1E6, t_timeDelay = -1E6;
    double t_elec_xs=0., t_elec_ys=0., t_elec_zs=0., t_elec_xe=0., t_elec_ye=0., t_elec_ze=0.;
    double t_und_xs=0., t_und_ys=0., t_und_zs=0., t_und_xe=0., t_und_ye=0., t_und_ze=0.;
    double t_elec_energy = 0., t_elec_time = 0., t_und_time = 0.;
   	int t_nb_cell_elec = 0, t_nb_cell_alpha = 0;
    bool t_elec_crossing;
    double t_eID, t_undID ;
    double t_nb_elec_crossing;
    std::vector<vector<double>> t_distance_elec_UND_per_elec_crossing , t_delta_t_cells_of_UND_per_elec_crossing, t_delta_t_cells_of_elec_per_elec_crossing;
    int t_nb_of_elec_candidates;
    int t_nPairs = 0;
    double t_mean_delta_t_UND_per_elec_crossing = std::numeric_limits<double>::quiet_NaN();
    double time_of_the_run = 0.0;
    double run_start_time = 0.0;
    int t_elec_is_associated_with_track;
    int run_number = run;
    int OM_num_elec; 

    double x_reference_source_elec;
    double y_reference_source_elec;
    double z_reference_source_elec;
    double x_reference_source_und;
    double y_reference_source_und;
    double z_reference_source_und;
    double vertex_is_on_source_foil_elec;
    double vertex_is_on_source_foil_und;
    
    double xM, yM, zM;
    double distance_vertex;


    double distance_gamma_BiPo_vertex = -1; 
    vector<double> delta_t_e_gamma_th_vs_meas;
    vector<double> E_gamma_BiPo; 
    int n_gamm_with_BiPo = 0;


    double c = 299792458; // speed of light in m/s
    double me = 9.1093837139E-31; // electron mass in kg

    pairTree->Branch("event", &t_event, "event/I");
    pairTree->Branch("elecIdx", &t_elecIdx, "elecIdx/I");
    pairTree->Branch("elecFit", &t_elecFit, "elecFit/I");
    pairTree->Branch("undIdx", &t_undIdx, "undIdx/I");
    pairTree->Branch("undFit", &t_undFit, "undFit/I");
    pairTree->Branch("minDist", &t_minDist, "minDist/D");
    pairTree->Branch("timeDelay", &t_timeDelay, "timeDelay/D");
    pairTree->Branch("elec_xs", &t_elec_xs, "elec_xs/D");
    pairTree->Branch("elec_ys", &t_elec_ys, "elec_ys/D");
    pairTree->Branch("elec_zs", &t_elec_zs, "elec_zs/D");
    pairTree->Branch("elec_xe", &t_elec_xe, "elec_xe/D");
    pairTree->Branch("elec_ye", &t_elec_ye, "elec_ye/D");
    pairTree->Branch("elec_ze", &t_elec_ze, "elec_ze/D");
    pairTree->Branch("und_xs", &t_und_xs, "und_xs/D");
    pairTree->Branch("und_ys", &t_und_ys, "und_ys/D");
    pairTree->Branch("und_zs", &t_und_zs, "und_zs/D");
    pairTree->Branch("und_xe", &t_und_xe, "und_xe/D");
    pairTree->Branch("und_ye", &t_und_ye, "und_ye/D");
    pairTree->Branch("und_ze", &t_und_ze, "und_ze/D");
    pairTree->Branch("elec_energy", &t_elec_energy, "elec_energy/D");
    pairTree->Branch("elec_time", &t_elec_time, "elec_time/D");
    pairTree->Branch("und_time", &t_und_time, "und_time/D");
	pairTree->Branch("nb_cell_elec", &t_nb_cell_elec, "nb_cell_elec/I");
    pairTree->Branch("nb_cell_alpha", &t_nb_cell_alpha);
    pairTree->Branch("elec_crossing", &t_elec_crossing);
    pairTree->Branch("elecClusterID", &t_elecClusterID, "elecClusterID/I");
    pairTree->Branch("undClusterID", &t_undClusterID, "undClusterID/I");
    pairTree->Branch("elec_is_associated_with_track", &t_elec_is_associated_with_track);
    pairTree->Branch("OM_num_elec", &OM_num_elec);

    pairTree->Branch("x_reference_source_elec", &x_reference_source_elec);
    pairTree->Branch("y_reference_source_elec", &y_reference_source_elec);
    pairTree->Branch("z_reference_source_elec", &z_reference_source_elec);

    pairTree->Branch("x_reference_source_und", &x_reference_source_und);
    pairTree->Branch("y_reference_source_und", &y_reference_source_und);
    pairTree->Branch("z_reference_source_und", &z_reference_source_und);
    pairTree->Branch("vertex_is_on_source_foil_elec", &vertex_is_on_source_foil_elec);
    pairTree->Branch("vertex_is_on_source_foil_und", &vertex_is_on_source_foil_und);
    pairTree->Branch("xM", &xM);
    pairTree->Branch("yM", &yM);
    pairTree->Branch("zM", &zM);
    pairTree->Branch("distance_vertex", &distance_vertex);

    pairTree->Branch("delta_t_e_gamma_th_vs_meas", &delta_t_e_gamma_th_vs_meas);
    pairTree->Branch("E_gamma_BiPo", &E_gamma_BiPo);
    pairTree->Branch("n_gamm_with_BiPo", &n_gamm_with_BiPo);


    pairTree->Branch("distance_elec_UND_per_elec_crossing", &t_distance_elec_UND_per_elec_crossing);
    pairTree->Branch("nb_elec_crossing", &t_nb_elec_crossing);
    pairTree->Branch("nb_of_elec_candidates", &t_nb_of_elec_candidates);
    pairTree->Branch("nPairs", &t_nPairs, "nPairs/I");
    pairTree->Branch("mean_delta_t_UND_per_elec_crossing", &t_mean_delta_t_UND_per_elec_crossing);
    pairTree->Branch("delta_t_cells_of_UND_per_elec_crossing", &t_delta_t_cells_of_UND_per_elec_crossing);
    pairTree->Branch("delta_t_cells_of_elec_per_elec_crossing", &t_delta_t_cells_of_elec_per_elec_crossing);
    pairTree->Branch("time_of_the_run", &time_of_the_run);
    pairTree->Branch("run_start_time", &run_start_time);
    pairTree->Branch("run_number", &run_number);
    if(info.found){
        std::cout<<"Run found in database, start time (unix) = "<<info.run_start<<", duration = "<<info.duration<<" ns, trip time = "<<info.tripTime<<" ns"<<std::endl;
    }
    else{
        std::cout<<"Run not found in database, all events will be processed"<<std::endl;
    }

        // Loop over all entries in the tree.

    for(int entry = 0; entry <= particules->GetEntries(); entry++) {



        //sndisplay::demonstrator *sndemonstrator = new sndisplay::demonstrator ("demonstrator_test");
        //sndemonstrator->setrange(0,1);
        particules->GetEntry(entry);
        
        delta_t_e_gamma_th_vs_meas.clear();
        E_gamma_BiPo.clear();
        n_gamm_with_BiPo = 0;

        if(info.found){
            
            
            double time_entry;
            double stop_time = info.tripTime;
            double unix_time_to_stop = info.run_start + stop_time;
            
            
            if(OM_timestamp_per_elec_cluster->size()>0){
                time_entry= OM_timestamp_per_elec_cluster->at(0).at(0);
            }
            //std::cout<<"Event "<<event_number<<" time: "<<time_entry<<" & tunix_time_to_stop: "<<unix_time_to_stop<<std::endl;
            if(time_entry>= stop_time){
                
                std::cout<<"Event "<<event_number<<" skipped because after the stop time of the run ("<<stop_time<<")"<<std::endl;
                time_of_the_run = stop_time;
                
                break;  
            }
            time_of_the_run = stop_time;
        }
        else{
            time_of_the_run = info.duration;
        }
        run_start_time = info.run_start;

        

		if(tracks_with_associated_calo && tracks_without_associated_calo){
            //std::cout<<"EVENT NUMBER: "<<event_number<<std::endl;
		    if(nb_of_elec_candidates >=1 && nb_of_UND_candidates>=1){
		        
		        struct FitInfo {
		            int idx;       // candidate index (electron & UND)
		            int fitIdx;    // fit solution index
		            // vertex start / end  (use as needed)
		            double xs, ys, zs;
		            double xe, ye, ze;
		            double time = 0.0;
                    int cluster_ID;
                    int elec_track_with_cluster;
                    int OM_num;
                    int vertex_is_on_source_foil;
                    // Vertex position on reference source plane
                    // Initialize to NaN so we don't propagate uninitialized memory
                    double x_reference_source = std::numeric_limits<double>::quiet_NaN();
                    double y_reference_source = std::numeric_limits<double>::quiet_NaN();
                    double z_reference_source = std::numeric_limits<double>::quiet_NaN();

                    double und_elec_crossing_xs, und_elec_crossing_ys, und_elec_crossing_zs, und_elec_crossing_xe, und_elec_crossing_ye, und_elec_crossing_ze;

		        };
		        struct PairInfo {
		            int elecIdx, elecFit;
		            int undIdx, undFit;
		            double minDist;
		            double timeDelay;
                    int cluster_ID_elec;
                    int cluster_ID_und;
                    int OM_num;
                    // initialize to NaN by default to avoid garbage if not set
                    double x_reference_source_elec = std::numeric_limits<double>::quiet_NaN();
                    double y_reference_source_elec = std::numeric_limits<double>::quiet_NaN();
                    double z_reference_source_elec = std::numeric_limits<double>::quiet_NaN();
                    double x_reference_source_und = std::numeric_limits<double>::quiet_NaN();
                    double y_reference_source_und = std::numeric_limits<double>::quiet_NaN();
                    double z_reference_source_und = std::numeric_limits<double>::quiet_NaN();
                    int vertex_is_on_source_foil_elec;
                    int vertex_is_on_source_foil_und;
		        };
			
		        // collect fits
		        std::vector<FitInfo> elecFits;
		        std::vector<FitInfo> undFits;
			
		        for(int elec_idx = 0; elec_idx < nb_of_elec_candidates; ++elec_idx){
		            int n_elec_fit = nb_fit_solution_per_elec_cluster->at(elec_idx);
		            for(int elec_fit_idx = 0; elec_fit_idx < n_elec_fit; ++elec_fit_idx){
		                FitInfo f;
		                f.idx = elec_idx; f.fitIdx = elec_fit_idx;
		                f.xs = x_start_per_elec_cluster->at(elec_idx).at(elec_fit_idx);
		                f.ys = y_start_per_elec_cluster->at(elec_idx).at(elec_fit_idx);
		                f.zs = z_start_per_elec_cluster->at(elec_idx).at(elec_fit_idx);
		                f.xe = x_end_per_elec_cluster->at(elec_idx).at(elec_fit_idx);
		                f.ye = y_end_per_elec_cluster->at(elec_idx).at(elec_fit_idx);
		                f.ze = z_end_per_elec_cluster->at(elec_idx).at(elec_fit_idx);
		                f.time = OM_timestamp_per_elec_cluster->at(elec_idx).at(0);
                        f.cluster_ID = ID_clsuter_per_elec_cluster->at(elec_idx);
                        f.elec_track_with_cluster = elec_is_associated_with_track->at(elec_idx);
                        f.OM_num = OM_num_per_elec_cluster->at(elec_idx).at(0);
                        
                        f.x_reference_source = x_is_on_reference_source_plane->at(elec_idx).at(elec_fit_idx);
                        f.y_reference_source = y_is_on_reference_source_plane->at(elec_idx).at(elec_fit_idx);
                        f.z_reference_source = z_is_on_reference_source_plane->at(elec_idx).at(elec_fit_idx);
                        f.vertex_is_on_source_foil = vertex_is_on_source_foil_per_elec_cluster->at(elec_idx).at(elec_fit_idx);
                       
                        

                        
                        elecFits.push_back(f);
		            }
		        }
			
		        for(int und_idx = 0; und_idx < nb_of_UND_candidates; ++und_idx){
		            int n_und_fit = nb_fit_solution_per_UND_cluster->at(und_idx);
		            for(int und_fit_idx = 0; und_fit_idx < n_und_fit; ++und_fit_idx){
		                FitInfo u;
		                u.idx = und_idx; u.fitIdx = und_fit_idx;
		                u.xs = x_start_per_UND_cluster->at(und_idx).at(und_fit_idx);
		                u.ys = y_start_per_UND_cluster->at(und_idx).at(und_fit_idx);
		                u.zs = z_start_per_UND_cluster->at(und_idx).at(und_fit_idx);
		                u.xe = x_end_per_UND_cluster->at(und_idx).at(und_fit_idx);
		                u.ye = y_end_per_UND_cluster->at(und_idx).at(und_fit_idx);
		                u.ze = z_end_per_UND_cluster->at(und_idx).at(und_fit_idx);
		                u.time = min_time_per_UND_cluster->at(und_idx); 
                        u.cluster_ID = ID_clsuter_UND->at(und_idx);
                        u.x_reference_source = x_is_on_reference_source_plane_per_UND_cluster->at(und_idx).at(und_fit_idx);
                        u.y_reference_source = y_is_on_reference_source_plane_per_UND_cluster->at(und_idx).at(und_fit_idx);
                        u.z_reference_source = z_is_on_reference_source_plane_per_UND_cluster->at(und_idx).at(und_fit_idx);
                        u.vertex_is_on_source_foil = vertex_is_on_source_foil_per_UND_cluster->at(und_idx).at(und_fit_idx);
		                
						undFits.push_back(u);
		            }
		        }
			
		        // compute all pair distances with delta t 
		        std::vector<PairInfo> allPairs;
		        for(const auto &ef : elecFits){
		            for(const auto &uf : undFits){
		                double min_dist = find_min_distance_between_tracks(
		                    ef.idx, uf.idx,
		                    x_start_per_elec_cluster, y_start_per_elec_cluster, z_start_per_elec_cluster,
		                    x_end_per_elec_cluster,   y_end_per_elec_cluster,   z_end_per_elec_cluster,
		                    x_start_per_UND_cluster,  y_start_per_UND_cluster,  z_start_per_UND_cluster,
		                    x_end_per_UND_cluster,    y_end_per_UND_cluster,    z_end_per_UND_cluster,
		                    ef.fitIdx, uf.fitIdx
		                );
		                PairInfo p;
		                p.elecIdx = ef.idx; p.elecFit = ef.fitIdx;
		                p.undIdx  = uf.idx; p.undFit  = uf.fitIdx;
		                p.minDist = min_dist;
		                p.timeDelay = uf.time - ef.time; 
                        p.cluster_ID_elec = ef.cluster_ID;
                        p.cluster_ID_und = uf.cluster_ID;
                        p.OM_num = ef.OM_num;
                        
                        p.x_reference_source_elec = ef.x_reference_source;
                        p.y_reference_source_elec = ef.y_reference_source;
                        p.z_reference_source_elec = ef.z_reference_source;

                        p.x_reference_source_und = uf.x_reference_source;
                        p.y_reference_source_und = uf.y_reference_source;
                        p.z_reference_source_und = uf.z_reference_source;

                        p.vertex_is_on_source_foil_elec = ef.vertex_is_on_source_foil;
                        p.vertex_is_on_source_foil_und = uf.vertex_is_on_source_foil;
		                
                        allPairs.push_back(p);
		            }
		        }
			
		        // print all pairs to check ! :D  
		        //for(const auto &p : allPairs){
		        //    std::cout << "Entry " << entry
		        //              << " | Elec " << p.elecIdx << " (fit " << p.elecFit << ")"
		        //              << " <-> UND " << p.undIdx << " (fit " << p.undFit << ")"
		        //              << " : dist = " << p.minDist << " mm"
		        //              << " | dt = " << p.timeDelay << std::endl;
		        //}
			
		        //  Here we gonna find best match per electron
		        std::sort(allPairs.begin(), allPairs.end(), [](const PairInfo &a, const PairInfo &b){
		            return a.minDist < b.minDist;
		        });
		        std::vector<bool> elecTaken(nb_of_elec_candidates, false);
		        std::vector<bool> undTaken(nb_of_UND_candidates, false);
		        std::vector<PairInfo> bestMatches;
		        for(const auto &p : allPairs){
		            if(!elecTaken[p.elecIdx] && !undTaken[p.undIdx]){
		                
		                bestMatches.push_back(p);
		                elecTaken[p.elecIdx] = true;
		                undTaken[p.undIdx] = true;
		            }
		        }

                // if there are multiple pairs for this event, log it (we keep one row per pair)
                if(bestMatches.size() > 1){
                   // std::cout << "Event " << event_number << " has " << bestMatches.size() << " pairs. Writing one row per pair." << std::endl;
                }

                // fill the pairTree for each best match
                for(const auto &p : bestMatches){
                    // record how many pairs were produced for this event
                    t_nPairs = static_cast<int>(bestMatches.size());
         	        t_event = event_number;
                    
        		    t_elecIdx = p.elecIdx;
        		    t_elecFit = p.elecFit;
        		    t_undIdx = p.undIdx;
        		    t_undFit = p.undFit;
        		    t_minDist = p.minDist;
        		    t_timeDelay = p.timeDelay;
        		    t_elecClusterID = p.cluster_ID_elec;
        		    t_undClusterID = p.cluster_ID_und;
                    
        		    
        		    t_elec_xs = x_start_per_elec_cluster->at(p.elecIdx).at(p.elecFit);

                    x_reference_source_elec = p.x_reference_source_elec;
                    y_reference_source_elec = p.y_reference_source_elec;
                    z_reference_source_elec = p.z_reference_source_elec;
                    x_reference_source_und  = p.x_reference_source_und;
                    y_reference_source_und  = p.y_reference_source_und;
                    z_reference_source_und  = p.z_reference_source_und;
                    vertex_is_on_source_foil_elec = p.vertex_is_on_source_foil_elec;
                    vertex_is_on_source_foil_und = p.vertex_is_on_source_foil_und;

                    t_elec_ys = y_start_per_elec_cluster->at(p.elecIdx).at(p.elecFit);
        		    t_elec_zs = z_start_per_elec_cluster->at(p.elecIdx).at(p.elecFit);
        		    t_elec_xe = x_end_per_elec_cluster->at(p.elecIdx).at(p.elecFit);
        		    t_elec_ye = y_end_per_elec_cluster->at(p.elecIdx).at(p.elecFit);
        		    t_elec_ze = z_end_per_elec_cluster->at(p.elecIdx).at(p.elecFit);
        		    t_elec_energy = E_OM_per_elec_cluster->at(p.elecIdx).at(0);
        		    t_elec_time = OM_timestamp_per_elec_cluster->at(p.elecIdx).at(0);
                    OM_num_elec = p.OM_num;
                    
                    //if(std::find(OM2kill.begin(), OM2kill.end(), OM_num_elec) != OM2kill.end()){ // OM Flag as bad ? if yes then pass this candidates
                    //    // skip this pair but continue filling other pairs for the same event
                    //    continue;
                    //}
                    

                    
        		    // fill UND info 
        		    t_und_xs = x_start_per_UND_cluster->at(p.undIdx).at(p.undFit);
        		    t_und_ys = y_start_per_UND_cluster->at(p.undIdx).at(p.undFit);
        		    t_und_zs = z_start_per_UND_cluster->at(p.undIdx).at(p.undFit);
        		    t_und_xe = x_end_per_UND_cluster->at(p.undIdx).at(p.undFit);
        		    t_und_ye = y_end_per_UND_cluster->at(p.undIdx).at(p.undFit);
        		    t_und_ze = z_end_per_UND_cluster->at(p.undIdx).at(p.undFit);
        		    t_und_time = min_time_per_UND_cluster->at(p.undIdx); 
					t_nb_cell_elec = nb_cell_per_elec_cluster->at(p.elecIdx);
					t_nb_cell_alpha = nb_cell_per_UND_cluster->at(p.undIdx);
                    t_elec_crossing = elec_used->at(p.elecIdx);
                    t_nb_of_elec_candidates = nb_of_elec_candidates;
                    // Clear vectors before filling
                    t_distance_elec_UND_per_elec_crossing.clear();
                    t_delta_t_cells_of_UND_per_elec_crossing.clear();
                    t_delta_t_cells_of_elec_per_elec_crossing.clear();
                    t_mean_delta_t_UND_per_elec_crossing = 1000000;
                    // Other information: 
                    t_nb_elec_crossing = nb_elec_crossing;
                    t_elec_is_associated_with_track = elec_is_associated_with_track->at(p.elecIdx);

                    if(elec_used->at(p.elecIdx)){
                        t_distance_elec_UND_per_elec_crossing.push_back(distance_elec_UND_per_elec_crossing->at(p.elecIdx));
                        
                        t_delta_t_cells_of_UND_per_elec_crossing.push_back(delta_t_cells_of_UND_per_elec_crossing->at(p.elecIdx));
                        t_delta_t_cells_of_elec_per_elec_crossing.push_back(delta_t_cells_of_elec_per_elec_crossing->at(p.elecIdx));

                        const auto &vals = delta_t_cells_of_UND_per_elec_crossing->at(p.elecIdx);
                        if(!vals.empty()){
                            double sum = std::accumulate(vals.begin(), vals.end(), 0.0);
                            t_mean_delta_t_UND_per_elec_crossing = sum / vals.size();
                        } else {
                            t_mean_delta_t_UND_per_elec_crossing = std::numeric_limits<double>::quiet_NaN();
                        }
                    }
                    else{
                        t_distance_elec_UND_per_elec_crossing.push_back({-1});
                        
                        t_delta_t_cells_of_UND_per_elec_crossing.push_back({-1});
                        t_delta_t_cells_of_elec_per_elec_crossing.push_back({-1});
                        t_mean_delta_t_UND_per_elec_crossing =-1;

                    }
                    

                    // Half distance between the two tracks calculation : https://www.youtube.com/watch?v=HC5YikQxwZA :) 
                    


                    if(t_elec_crossing == false){
                        Vec3 A  = {t_elec_xs, t_elec_ys, t_elec_zs};
                        Vec3 B  = {t_elec_xe, t_elec_ye, t_elec_ze};
                        Vec3 A2 = {t_und_xs, t_und_ys, t_und_zs};
                        Vec3 B2 = {t_und_xe, t_und_ye, t_und_ze};

                        Vec3 C, C2;
                        if (closestPointsBetweenLines(A, B, A2, B2, C, C2)) {
                            //cout << "Point C  sur D  : (" << C.x  << ", " << C.y  << ", " << C.z  << ")\n";
                            //cout << "Point C' sur D' : (" << C2.x << ", " << C2.y << ", " << C2.z << ")\n";
                            
                            // Distance minimale entre les droites
                            double distance = norm(C - C2);
                            //cout << "Distance minimale = " << distance << endl;
                            
                            // Point du milieu du segment CC'
                            Vec3 M = (C + C2) * 0.5;
                            //cout << "Milieu M du segment CC' : (" 
                            //     << M.x << ", " << M.y << ", " << M.z << ")\n";
                            xM = M.x;
                            yM = M.y;
                            zM = M.z;
                            distance_vertex = distance;
                        }
                    
                    } 
                    
                    // Here to modified but first need to put correct vector size of und_elec_crossing_ze ! to be able to got with the good idx first & last position ! 
                    else if(t_elec_crossing == true){
                        Vec3 A  = {und_elec_crossing_xs->at(t_elecIdx).at(0), und_elec_crossing_ys->at(t_elecIdx).at(0), und_elec_crossing_zs->at(t_elecIdx).at(0)};
                        Vec3 B  = {und_elec_crossing_xe->at(t_elecIdx).at(0), und_elec_crossing_ye->at(t_elecIdx).at(0), und_elec_crossing_ze->at(t_elecIdx).at(0)};
                        Vec3 A2 = {t_und_xs, t_und_ys, t_und_zs};
                        Vec3 B2 = {t_und_xe, t_und_ye, t_und_ze};

                        Vec3 C, C2;
                        if (closestPointsBetweenLines(A, B, A2, B2, C, C2)) {
                            //cout << "Point C  sur D  : (" << C.x  << ", " << C.y  << ", " << C.z  << ")\n";
                            //cout << "Point C' sur D' : (" << C2.x << ", " << C2.y << ", " << C2.z << ")\n";
                            
                            // Distance minimale entre les droites
                            double distance = norm(C - C2);
                            //cout << "Distance minimale = " << distance << endl;
                            
                            // Point du milieu du segment CC'
                            Vec3 M = (C + C2) * 0.5;
                            //cout << "Milieu M du segment CC' : (" 
                            //     << M.x << ", " << M.y << ", " << M.z << ")\n";
                            xM = M.x;
                            yM = M.y;
                            zM = M.z;
                            distance_vertex = distance;
                        }
                    }

                    double elec_track_length = -1;
                    // Compute electron track length
                    if(t_elec_crossing == false){
                        if(OM_num_elec<520){

                            elec_track_length = sqrt( pow((x_is_on_main_calorimeter->at(t_elecIdx).at(t_elecFit) - t_elec_xs),2) + pow((y_is_on_main_calorimeter->at(t_elecIdx).at(t_elecFit) - t_elec_ys),2) + pow((z_is_on_main_calorimeter->at(t_elecIdx).at(t_elecFit) - t_elec_zs),2) )/1000; // in meters
                            //std::cout<<"Electron track length: "<<elec_track_length<<std::endl;

                        }
                        else if(OM_num_elec>=520){
                            elec_track_length = sqrt( pow((x_is_on_x_calorimeter->at(t_elecIdx).at(t_elecFit) - t_elec_xs),2) + pow((y_is_on_x_calorimeter->at(t_elecIdx).at(t_elecFit) - t_elec_ys),2) + pow((z_is_on_x_calorimeter->at(t_elecIdx).at(t_elecFit) - t_elec_zs),2) )/1000; // in meters
                        }
                        
                    }
                    double second_part_elec_length;
                    
                    if(t_elec_crossing == true){
                        double firt_part_elec_length = sqrt( pow((x_reference_source_elec - t_elec_xs),2) + pow((y_reference_source_elec - t_elec_ys),2) + pow((z_reference_source_elec - t_elec_zs),2) ); // in mm
                        if(OM_num_elec<520){
                            second_part_elec_length = sqrt( pow((x_is_on_main_calorimeter->at(t_elecIdx).at(t_elecFit) - x_reference_source_elec),2) + pow((y_is_on_main_calorimeter->at(t_elecIdx).at(t_elecFit) - y_reference_source_elec),2) + pow((z_is_on_main_calorimeter->at(t_elecIdx).at(t_elecFit) - z_reference_source_elec),2) ); // in mm    
                        }
                        else if(OM_num_elec>=520){
                            second_part_elec_length = sqrt( pow((x_is_on_x_calorimeter->at(t_elecIdx).at(t_elecFit) - x_reference_source_elec),2) + pow((y_is_on_x_calorimeter->at(t_elecIdx).at(t_elecFit) - y_reference_source_elec),2) + pow((z_is_on_x_calorimeter->at(t_elecIdx).at(t_elecFit) - z_reference_source_elec),2) ); // in mm
                        }
                        

                        elec_track_length = (firt_part_elec_length + second_part_elec_length)/1000; // in meters
                    }

                    if(nb_isolated_calo > 0){
                        for(int calo_idx =0; calo_idx < nb_isolated_calo; calo_idx++){
                            double calo_x = x_calo_gamma->at(calo_idx);
                            double calo_y = y_calo_gamma->at(calo_idx);
                            double calo_z = z_calo_gamma->at(calo_idx);
                            distance_gamma_BiPo_vertex = -1; 
                            if(std::isnan(xM) || std::isnan(yM) || std::isnan(zM)){
                                distance_gamma_BiPo_vertex = sqrt( pow((calo_x - t_elec_xs),2) + pow((calo_y - t_elec_ys),2) + pow((calo_z - t_elec_zs),2) )/1000; // in meters
                            }
                            else{
                                distance_gamma_BiPo_vertex = sqrt( pow((calo_x - xM),2) + pow((calo_y - yM),2) + pow((calo_z - zM),2) )/1000; // in meters
                                //std::cout<<"Distance calo to M: "<<distance_calo_to_M<<std::endl;
                            }
                            double expected_TOF_gamma = distance_gamma_BiPo_vertex/c;
                            double measured_TOF_gamma = t_elec_time - isolated_calo_timestamp->at(calo_idx);
                            double beta = sqrt(t_elec_energy*(t_elec_energy + 2*0.511/pow(c, 2)))/(t_elec_energy+ 2*0.511);
                            double expected_TOF_elec = elec_track_length/(beta*c);
                            n_gamm_with_BiPo ++;

                            //delta_t_e_gamma_th_vs_meas.push_back(expected_TOF_gamma - measured_TOF_gamma);
                            delta_t_e_gamma_th_vs_meas.push_back(isolated_calo_timestamp->at(calo_idx) - t_elec_time - (expected_TOF_gamma - expected_TOF_elec) );
                            E_gamma_BiPo.push_back(E_isolated_calo->at(calo_idx));
                        }
                    }
                    
                    pairTree->Fill();
        		}
		    }
		}		
	}
    outFile->cd();
	outFile->Write();
	outFile->Close();
	std::cout << "Root file saved !" << std::endl;
         
    //std::string path = "/home/antoine/Documents/These/1-Radon/PTD/data_extracted_BiPo_or_not/img_PTD_extracted/" + std::to_string(event_number) + ".png";
    //sndemonstrator->canvas->SaveAs(path.c_str());
    //sndemonstrator->canvas->Close();
      

}
