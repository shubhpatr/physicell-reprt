/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = basic_volume_model; // standard_volume_update_function;
//	cell_defaults.functions.update_velocity = standard_update_cell_velocity;
	cell_defaults.functions.update_velocity = heterotypic_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = plasto_elastic_mechanics; // NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 

	// let's add custom data not in the XML 
	std::vector<double> zero = {0,0,0};
	cell_defaults.custom_data.add_vector_variable( "rest_position" , "microns" , zero ); 
	cell_defaults.custom_data.add_vector_variable( "BM_attach_point" , "microns", zero ); 

	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 
	
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_subcells_from_pugixml(); 	

	// more custom setup 
	// // for each cell, set its relaxed position to current position 
	int nRP = 0; // cell_defaults.custom_data.find_vector_variable( "rest_position");
	for( int n=0 ; n < (*all_cells).size(); n++ )
	{
		Cell* pC = (*all_cells)[n]; 
		pC->custom_data.vector_variables[nRP].value = pC->position; 
	}
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring_wrapped(pCell); } 
// paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 


bool load_subcells_from_pugixml( pugi::xml_node root )
{
	pugi::xml_node node = root.child( "initial_conditions" ); 
	if( !node )
	{ 
		std::cout << "Warning: XML-based cell positions has wrong formating. Ignoring!" << std::endl; 
		return false;
	}

	node = node.child( "cell_positions" ); 
	if( !node )
	{
		std::cout << "Warning: XML-based cell positions has wrong formating. Ignoring!" << std::endl; 
		 return false;
	}

	// enabled? 
	if( node.attribute("enabled").as_bool() == false )
	{ return false; }

	// get filename 

	std::string folder = xml_get_string_value( node, "folder" ); 
	std::string filename = xml_get_string_value( node, "filename" ); 
	std::string input_filename = folder + "/" + filename; 

	std::string filetype = node.attribute("type").value() ; 

	// what kind? 
	if( filetype == "csv" || filetype == "CSV" )
	{
		std::cout << "Loading cells from CSV file " << input_filename << " ... " << std::endl; 
		load_subcells_csv( input_filename );
		system("timeout 1");
		return true; 
	}
	if( filetype == "matlab" || filetype == "mat" || filetype == "MAT" )
	{
		std::cout << "Error: Load cell positions from matlab not yet supported. Try CSV." << std::endl; 
		exit(-1); 
		std::cout << "Loading cells from matlab file " << input_filename << " ... " << std::endl; 
		return false; 
	}
	if( filetype == "scene" )
	{
		std::cout << "Error: load cell positions from scene not yet supported. Try CSV." << std::endl; 
		exit(-1); 
		std::cout << "Loading cells from scene file " << input_filename << " ... " << std::endl; 
		return false; 
	}
	if( filetype == "physicell" || filetype == "PhysiCell" )
	{
		std::cout << "Error: load cell positions from PhysiCell snapshot not yet supported. Try CSV." << std::endl; 
		exit(-1); 
		std::cout << "Loading cells from PhysiCell file " << input_filename << " ... " << std::endl; 
		return false; 
	}

	return false; 
}

void load_subcells_csv( std::string filename )
{
	std::ifstream file( filename, std::ios::in );
	if( !file )
	{ 
		std::cout << "Error: " << filename << " not found during cell loading. Quitting." << std::endl; 
		exit(-1);
	}
	
	static int nCellID = cell_defaults.custom_data.find_variable_index( "cell_ID" ); 

	std::string line;
	while (std::getline(file, line))
	{
		std::vector<double> data;
		csv_to_vector( line.c_str() , data ); 

		if( data.size() != 5 )
		{
			std::cout << "Error! Importing subcells from a CSV file expects each row to be x,y,z,typeID,cellID." << std::endl;
			exit(-1);
		}

		std::vector<double> position = { data[0] , data[1] , data[2] };

		int my_type = (int) data[3]; 
		int my_ID = (int) data[4]; 
		Cell_Definition* pCD = find_cell_definition( my_type );
		if( pCD != NULL )
		{
			std::cout << "Creating " << pCD->name << " (type=" << pCD->type << ") at " << position << std::endl; 
			Cell* pCell = create_cell( *pCD ); 
			pCell->assign_position( position ); 
			
			pCell->custom_data[nCellID] = my_ID; 
		}
		else
		{
			std::cout << "Warning! No cell definition found for index " << my_type << "!" << std::endl
			<< "\tIgnoring cell in " << filename << " at position " << position << std::endl; 
		}

	}

	file.close(); 	
}

bool load_subcells_from_pugixml( void )
{ return load_subcells_from_pugixml( physicell_config_root ); }

std::vector<std::string> paint_by_number_cell_coloring_wrapped( Cell* pCell )
{
	static std::vector< std::string > colors(0); 
	static bool setup_done = false; 
	if( setup_done == false )
	{
		colors.push_back( "grey" ); // default color will be grey 

		colors.push_back( "red" );
		colors.push_back( "yellow" ); 	
		colors.push_back( "green" ); 	
		colors.push_back( "blue" ); // 4
		
		colors.push_back( "magenta" ); 	
		colors.push_back( "orange" ); 	
		colors.push_back( "lime" ); 	
		colors.push_back( "cyan" ); // 8
		
		colors.push_back( "hotpink" ); 	
		colors.push_back( "peachpuff" ); 	
		colors.push_back( "darkseagreen" ); 	
		colors.push_back( "lightskyblue" ); // 12

		colors.push_back( "darkred" );
		colors.push_back( "goldenrod" ); 	
		colors.push_back( "darkgreen" ); 	
		colors.push_back( "darkblue" );  // 16
		
		setup_done = true; 
	}
	
	int number_of_colors = colors.size(); 
	
	// start all black 
	
	std::vector<std::string> output = { "black", "black", "black", "black" }; 
	
	// paint by number -- by cell type 
	
	// modular arithmetic  
	// int type = (int) ( ((unsigned int) pCell->type) % 13 ); 
	static int nCellID = cell_defaults.custom_data.find_variable_index( "cell_ID" ); 
	// int agent_type = pCell->type % number_of_colors; 
	int agent_type = (int) (pCell->custom_data[nCellID]) % number_of_colors; 
	
//	std::cout << type << std::endl; 
	
	std::string interior_color = "white"; 
	if( agent_type < number_of_colors )
	{ interior_color = colors[ agent_type ]; }
	
	output[0] = interior_color; // set cytoplasm color 
	
	if( pCell->phenotype.death.dead == false ) // if live, color nucleus same color 
	{
		output[2] = interior_color; 
		output[3] = interior_color; 
	}
	else
	{
		// apoptotic cells will retain a black nucleus 
		// if necrotic, color the nucleus brown 
		if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
		{
			output[2] = "rgb(139,69,19)";
			output[3] = "rgb(139,69,19)";
		}
	}
	
	return output; 
}
