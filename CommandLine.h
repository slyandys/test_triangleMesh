#ifndef _COMMAND_LINE_H_
#define _COMMAND_LINE_H_

#include <string>
#include <sstream>

using std::string;
using std::stringstream;

struct CommandLine {
  enum Variant { PARAM_ESTIMATE, TRACK, FEM_TRACK };
  Variant d_variant;
  string d_path; // template file -- called path for backwards compatibility
  string d_filenames;
  bool d_useVrml;
  string d_outputName;
  string d_contactFileName;
  float d_sizeNbhd;
  string d_lowResTetFile;
  int d_numMarkers;
  string d_markerFile;
  string d_forceFileName;
  string d_supportSurfFileName;
  string d_solidmeshFileName;
  string d_probeFileName;
  string d_anchorPointsFile;

  inline CommandLine( int argc, char *argv[] );

  inline void printHelp( std::ostream& os );
};


CommandLine::CommandLine( int argc, char *argv[] ) : d_variant( FEM_TRACK ),
  d_path(""), d_filenames(""), d_useVrml(false), d_outputName(""), d_contactFileName(""),
  d_sizeNbhd(-1), d_lowResTetFile(""), d_numMarkers(-1), d_markerFile(""), d_forceFileName(""),
  d_supportSurfFileName(""), d_solidmeshFileName(""), d_probeFileName(""), d_anchorPointsFile("")
{
  stringstream sstr;
  string token;
  for ( int i=1; i<argc; ++i) {
    sstr << argv[i] << ' ';
  }
  // get tokens
  while ( sstr >> token ) {
    if (  token[0] == '-' ) {
      // cerr << "Token: " << token << endl;
      // new parameter
      if ( token.compare( "-paramestimate" ) == 0 ) {
	d_variant = PARAM_ESTIMATE;
	continue;
      }
      if ( token.compare( "-track" ) == 0 ) {
	d_variant = TRACK;
	continue;
      }
      if ( token.compare( "-femtrack" ) == 0 ) {
	d_variant = FEM_TRACK;
	continue;
      }
      if ( token.compare( "-t" ) == 0 ) {
	sstr >> d_path;
	continue;
      }
      if ( token.compare( "-fn" ) == 0 ) {
	sstr >> d_filenames;
	continue;
      }
      if ( token.compare( "-vrml" ) == 0 ) {
	d_useVrml = true;
	continue;
      }
      if ( token.compare( "-o" ) == 0 ) {
	sstr >> d_outputName;
	continue;
      }
      if ( token.compare( "-c" ) == 0 ) {
	sstr >> d_contactFileName;
	continue;
      }
      if ( token.compare( "-ff" ) == 0 ) {
	sstr >> d_forceFileName;
	continue;
      }
      if ( token.compare( "-sf" ) == 0 ) {
	sstr >> d_supportSurfFileName;
	continue;
      }
      if ( token.compare( "-p" ) == 0 ) {
	sstr >> d_probeFileName;
	continue;
      }
      if ( token.compare( "-szn" ) == 0 ) {
	continue;
      }
      if ( token.compare( "-l" ) == 0 ) {
	sstr >> d_lowResTetFile;
	continue;
      }
      if ( token.compare( "-m" ) == 0 ) {
	sstr >> d_markerFile;
	sstr >> d_numMarkers;
	continue;
      }
      if ( token.compare( "-s" ) == 0 ) {
	sstr >> d_solidmeshFileName;
	continue;
      }
      if ( token.compare( "-a" ) == 0 ) {
	sstr >> d_anchorPointsFile;
	continue;
      }
      if ( token.compare( "-h" ) == 0 ) {
	printHelp( std::cerr );
	exit(0);
      }
    } else {
      std::cerr << "Unexpected parameter: " << token << std::endl;
      printHelp( std::cerr );
    }
  }
}

void CommandLine::printHelp( std::ostream& os ) {
  os << "tracker " << "[-paramestimate|-track|-femtrack] " << "-t templateFile " << "-fn filenames " <<
	  "{-vrml} " << "-o outputFile " << "[[-c contactFileName {-ff forceFileName} {-ss d_supportSurfFileName}]|-p probeFileName] " <<
    "-szn 3 " << "-l lowResTetFile " << "-m markerFile numMarkers " << "-s solidmeshFileName " <<
    "-a anchorPointsFile" << endl;
  return;
}





#endif





