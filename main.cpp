/* 
 * File:   main.cpp
 * Author: tombieli
 *
 * Created on 1 August 2016, 15:05
 */

#include <iostream>
#include <sstream>
#include <cstring>

#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <boost/program_options.hpp>
#include "correction.h"

using namespace std;

boost::program_options::variables_map init(int argc, char ** argv)
{
    boost::program_options::options_description options("GEOS Height correction options");
    options.add_options()
    ("help","Shows help message")
    ("input", boost::program_options::value<std::string>(), "Location of image containing height information.")
    ("output", boost::program_options::value<std::string>(), "Location of result raster with corrected coordiantes")
    ("height-band",boost::program_options::value<int>()->default_value(1),"Band with height information [m]. ")
    ("requierd-accuracy",boost::program_options::value<double>()->default_value(10),"Required accuracy [m].")
    ("iterations-limit",boost::program_options::value<int>()->default_value(100),"Maximum number of iteration of netwon method per pixel."); 

    auto ret = boost::program_options::variables_map();

    auto style = boost::program_options::command_line_style::unix_style;
    style = (decltype(style))(style ^ boost::program_options::command_line_style::allow_short);
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options, style), ret);
    boost::program_options::notify(ret);

    if(argc <= 1 || ret.count("help") >=1)
    {
        std::cout << options << std::endl;
        exit(1);
    }

    return ret;
}

void CopyGeoMetadata(GDALDataset* dest, GDALDataset* source)
{
    double geo[6];
    source->GetGeoTransform(geo);
    dest->SetGeoTransform(geo);
    dest->SetProjection(source->GetProjectionRef());
}



int main(int argc, char** argv) {
    auto variablesMap = init(argc,argv);
    GDALAllRegister();
    
    int ret = 0;
    
    GDALDataset *input = nullptr, *output = nullptr;
    
    try
    {
        auto inputPath = variablesMap["input"].as<std::string>();

        input = (GDALDataset*)GDALOpen(inputPath.c_str(), GA_ReadOnly);
        if(input == nullptr)
            throw runtime_error("Cannot open "+inputPath);
        
        auto tifDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
        if(tifDriver == nullptr)
            throw runtime_error("There is no GTiff driver");
        
        GDALDataType type = GDT_Float64;
        OGRSpatialReference srs;
        auto srs_c = input->GetProjectionRef();
        srs.importFromWkt(&srs_c);
        
        auto projectionName = string(srs.GetAttrValue("PROJECTION"));
        
        if(projectionName != "Geostationary_Satellite")
            throw runtime_error("Height correction requires Geostationary_Satellite projection.");
        

        auto outputPath = variablesMap["output"].as<std::string>();   
        
        output = tifDriver->Create(outputPath.c_str(),input->GetRasterXSize(), input->GetRasterYSize(),2, type, nullptr);
        
        CopyGeoMetadata(output,input);
        
        
        performCorrection(input,output,variablesMap["height-band"].as<int>(),srs,variablesMap["requierd-accuracy"].as<double>(),variablesMap["iterations-limit"].as<int>());
        
    }
    catch(exception &ex)
    {
        cerr << "Error: "<< ex.what() << endl;
        ret = -1;
    }
    
    if(input != nullptr)
    {
        GDALClose(input);
    }
    
    if(output != nullptr)
    {
        GDALClose(output);
    }
    return 0;
}

