#include <cassert>
#include <string>
#include <format>
#include "mcnptools/Ptrac.hpp"

#include "mcnptools/StringOps.hpp"

namespace mcnptools {

Ptrac::Ptrac(const std::string& filename, const unsigned int format):
  m_filename( filename ),
  m_format( static_cast<Ptrac::PtracFormat>(format) ) {

  if( format == Ptrac::HDF5_PTRAC ) {
    m_handle = {};
    auto file_h5 = f5::File(m_filename, 'r');
    m_hdf5_parser = std::make_unique<Ptrac::HDF5Parser>( MakeHDF5PtrackParser(file_h5, "ptrack") );
  }
  else
  {
    if( format == Ptrac::BIN_PTRAC ) {
      m_handle.open(filename.c_str(), std::ifstream::binary);
      if (m_handle.fail()) {
        throw McnpToolsException( std::format("Failed to open binary PTRAC file {}", filename));
      }
    }
    else { //ASCII format 
      m_handle.open(filename.c_str());
      if (m_handle.fail()) {
        throw McnpToolsException( std::format("Failed to open ASCII PTRAC file {}",filename) );
      }
    }
    ReadHeader();
  }

  // --- NEW: Pre-calculate Event Data Layouts ---
    auto merge_layout = [&](int event_type, LineIndex idx1, LineIndex idx2) {
        std::vector<int> layout;
        layout.reserve(m_datent[idx1].size() + m_datent[idx2].size());
        layout.insert(layout.end(), m_datent[idx1].begin(), m_datent[idx1].end());
        layout.insert(layout.end(), m_datent[idx2].begin(), m_datent[idx2].end());
        m_event_layouts[event_type] = std::move(layout);
    };

    merge_layout(Ptrac::SRC, IDX_SRC1, IDX_SRC2);
    merge_layout(Ptrac::BNK, IDX_BNK1, IDX_BNK2);
    merge_layout(Ptrac::SUR, IDX_SUR1, IDX_SUR2);
    merge_layout(Ptrac::COL, IDX_COL1, IDX_COL2);
    merge_layout(Ptrac::TER, IDX_TER1, IDX_TER2);
}

constexpr Ptrac::LineIndex GetBaseIndex(const std::string& type) {
    if (type == "src") return Ptrac::IDX_SRC1;
    if (type == "bnk") return Ptrac::IDX_BNK1;
    if (type == "sur") return Ptrac::IDX_SUR1;
    if (type == "col") return Ptrac::IDX_COL1;
    if (type == "ter") return Ptrac::IDX_TER1;
    return Ptrac::IDX_NPS;
}

void Ptrac::ReadHeader() {
  if( m_format == Ptrac::BIN_PTRAC ) {
    int size1, size2;

    // determine file size
    m_handle.seekg(0,std::ios::end);
    int64_t fsize = m_handle.tellg();
    m_handle.seekg(0,std::ios::beg);

    // read the version
    m_handle.read( (char*) &size1, sizeof(int) );

    if( size1 >= fsize || size1 != sizeof(int) ) {
      throw McnpToolsException( "Failed to read binary PTRAC" );
    }

    m_handle.read( (char*) &m_version, size1);
    m_handle.read( (char*) &size2, sizeof(int) );

    if( size1 != size2 || m_version != -1 ) {
      throw McnpToolsException( "Failed to read binary PTRAC" );
    }

    // read the code data 
    m_handle.read( (char*) &size1, sizeof(int) );
    
    char code[8], ver[5], loddat[8], idtm[19];
    m_handle.read( code, sizeof(code) );
    m_code = std::string(code,sizeof(code));
    m_handle.read( ver, sizeof(ver) );
    m_codever = std::string(ver,sizeof(ver));
    m_handle.read( loddat, sizeof(loddat) );
    m_loddat = std::string(loddat,sizeof(loddat));
    m_handle.read( idtm, sizeof(idtm) );
    m_idtm = std::string(idtm,sizeof(idtm));
    cjsoft::stringops::trim(m_idtm);

    m_handle.read( (char*) &size2, sizeof(int) );

    if( size1 != size2 ) {
      std::stringstream ss;
      ss << "Failed to read binary PTRAC";
      throw McnpToolsException( ss.str() );
    }

    // read the comment line
    m_handle.read( (char*) &size1, sizeof(int) );
    assert( size1 == 80 || size1 == 128 );

    m_comment.resize(size1);
    m_handle.read(&m_comment[0], size1*sizeof(char));
    
    m_handle.read( (char*) &size2, sizeof(int) );

    if( size1 != size2 ) {
      throw McnpToolsException( "Failed to read binary PTRAC" );
    }

    // read the keyword entries
    bool done = false;
    std::vector<double> kwent;
    unsigned int nkw = 0;
    unsigned int kwlinecnt = 0;
    while (! done) {
      kwlinecnt++;
      m_handle.read( (char*) &size1, sizeof(int) );

      double buffer[10];
      m_handle.read( (char*) buffer, 10*sizeof(double) );
      if( kwlinecnt == 1 ) {
        nkw = (unsigned int) buffer[0];
        kwent.insert(kwent.end(), &buffer[1], &buffer[1] + 9);
      }
      else {
        kwent.insert(kwent.end(), buffer, buffer + 10);
      }

      unsigned int nkwcnt = 0;
      unsigned int i=0;
      while( i<kwent.size() ) {
        nkwcnt++;
        unsigned int num_ent = (unsigned int) kwent[i];
        i += num_ent+1;
      }

      m_handle.read( (char*) &size2, sizeof(int) );

      if( size1 != size2 ) {
        throw McnpToolsException( "Failed to read binary PTRAC" );
      }

      if( nkwcnt >= nkw )
        done = true;
    }

    // read number data
    m_handle.read( (char*) &size1, sizeof(int) );

    int nnps, ipt, single_double, unused[7];
    int64_t nsrc1, nsrc2, nbnk1, nbnk2, nsur1, nsur2, ncol1, ncol2, nter1, nter2;
    m_handle.read( (char*) &nnps, sizeof(int) );
    m_handle.read( (char*) &nsrc1, sizeof(nsrc1) );
    m_handle.read( (char*) &nsrc2, sizeof(nsrc2) );
    m_handle.read( (char*) &nbnk1, sizeof(nbnk1) );
    m_handle.read( (char*) &nbnk2, sizeof(nbnk2) );
    m_handle.read( (char*) &nsur1, sizeof(nsur1) );
    m_handle.read( (char*) &nsur2, sizeof(nsur2) );
    m_handle.read( (char*) &ncol1, sizeof(ncol1) );
    m_handle.read( (char*) &ncol2, sizeof(ncol2) );
    m_handle.read( (char*) &nter1, sizeof(nter1) );
    m_handle.read( (char*) &nter2, sizeof(nter2) );
    m_handle.read( (char*) &ipt, sizeof(int) );
    m_handle.read( (char*) &single_double, sizeof(int) );
    m_handle.read( (char*) &unused, sizeof(unused) );

    m_nument[IDX_NPS] = nnps;
    m_nument[IDX_SRC1] = nsrc1;
    m_nument[IDX_SRC2] = nsrc2;
    m_nument[IDX_BNK1] = nbnk1;
    m_nument[IDX_BNK2] = nbnk2;
    m_nument[IDX_SUR1] = nsur1;
    m_nument[IDX_SUR2] = nsur2;
    m_nument[IDX_COL1] = ncol1;
    m_nument[IDX_COL2] = ncol2;
    m_nument[IDX_TER1] = nter1;
    m_nument[IDX_TER2] = nter2;


    m_handle.read( (char*) &size2, sizeof(int) );
    
    if( size1 != size2 ) {
      throw McnpToolsException( "Failed to read binary PTRAC" );
    }

    // read data types
    m_handle.read( (char*) &size1, sizeof(int) );

    for(auto x : m_lines) {
      for(unsigned int j=0; j<m_nument[ x ]; j++) {
        if( x == IDX_NPS ) {
          int64_t tmp;
          m_handle.read( (char*) &tmp, sizeof(tmp) );
          m_datent[ x ].push_back( tmp );
        }
        else {
          int tmp;
          m_handle.read( (char*) &tmp, sizeof(tmp) );
          m_datent[ x ].push_back( tmp );
        }
      }
    }

    m_handle.read( (char*) &size2, sizeof(int) );

    if( size1 != size2 ) {
      throw McnpToolsException( "Failed to read binary PTRAC" );
    }

  }
  else {
    // read the version
    m_handle >> m_version;

    // read the code data, if it is available (prdmp third entry effects this)
    std::string idtm1, idtm2;
    std::string optional_version_line;
    std::getline(m_handle, optional_version_line); // processes to the next line
    std::getline(m_handle, optional_version_line); 
    if (optional_version_line != " ") { 
      std::stringstream ss {optional_version_line};
      ss >> m_code >> m_codever >> m_loddat >> idtm1 >> idtm2;
      m_idtm = idtm1 + " " + idtm2;
    } // else, the above variables are default-constructed empty strings

    // read the comment line
    cjsoft::stringops::getline(m_handle, m_comment); // processes to the next line
    cjsoft::stringops::getline(m_handle, m_comment);

    // read the keyword entries
    bool done = false;
    std::vector<double> kwent;
    unsigned int nkw = 0;
    unsigned int kwlinecnt = 0;
    while (! done) {
      kwlinecnt++;

      if( kwlinecnt == 1 ) {
        double tmp;
        m_handle >> tmp;
        nkw = (unsigned int) tmp;
        for(unsigned int i=1; i<10; i++) {
          m_handle >> tmp;
          kwent.push_back(tmp);
        }
      }
      else {
        double tmp;
        for(unsigned int i=0; i<10; i++) {
          m_handle >> tmp;
          kwent.push_back(tmp);
        }
      }

      unsigned int nkwcnt = 0;
      unsigned int i=0;
      while( i<kwent.size() ) {
        nkwcnt++;
        unsigned int num_ent = (unsigned int) kwent[i];
        i += num_ent+1;
      }

      if( nkwcnt >= nkw )
        done = true;
    }

    // read number data
    int nnps, ipt, single_double, unused[7];
    int64_t nsrc1, nsrc2, nbnk1, nbnk2, nsur1, nsur2, ncol1, ncol2, nter1, nter2;

    m_handle >> nnps >> nsrc1 >> nsrc2 >> nbnk1 >> nbnk2 >> nsur1 >> nsur2 >> ncol1 >> ncol2 >> nter1 >> nter2 >> ipt >> single_double;
    for(unsigned int i=0; i<7; i++) {
      m_handle >> unused[i];
    }

    m_nument[IDX_NPS] = nnps;
    m_nument[IDX_SRC1] = nsrc1;
    m_nument[IDX_SRC2] = nsrc2;
    m_nument[IDX_BNK1] = nbnk1;
    m_nument[IDX_BNK2] = nbnk2;
    m_nument[IDX_SUR1] = nsur1;
    m_nument[IDX_SUR2] = nsur2;
    m_nument[IDX_COL1] = ncol1;
    m_nument[IDX_COL2] = ncol2;
    m_nument[IDX_TER1] = nter1;
    m_nument[IDX_TER2] = nter2;

    // read data types

    for(auto x : m_lines) {
      for(unsigned int j=0; j<m_nument[ x ]; j++) {
        if( x == IDX_NPS ) {
          int64_t tmp;
          m_handle >> tmp;
          m_datent[ IDX_NPS ].push_back( tmp );
        }
        else {
          int tmp;
          m_handle >> tmp;
          m_datent[ x ].push_back( tmp );
        }
      }
    }
  }
}

PtracHistory Ptrac::ReadHistory() {
  int size1, size2;

  PtracHistory hist;
  static_assert(std::is_nothrow_move_constructible<PtracHistory>::value);
  static_assert(std::is_nothrow_move_constructible<PtracEvent>::value);
  // read the nps line
  double next_event_type;
  if( m_format == Ptrac::BIN_PTRAC ) ReadValue(size1);

  PtracNps nps;
  // Reuse a buffer for binary reads to avoid repeated heap allocations
  std::vector<double> binary_buffer;
  for(unsigned int i=0; i<m_nument[IDX_NPS]; i++) {
    int64_t tmp;

    ReadValue(tmp);

    if( ! m_handle.good() )
      return hist;

    switch( m_datent[IDX_NPS][i] ) {
      case Ptrac::NPS:
        nps.m_nps = tmp;
        break;
      case Ptrac::FIRST_EVENT_TYPE:
        next_event_type = tmp;
        break;
      case Ptrac::NPSCELL:
        nps.m_cell = tmp;
        break;
      case Ptrac::NPSSURFACE:
        nps.m_surface = tmp;
        break;
      case Ptrac::TALLY:
        nps.m_tally = tmp;
        break;
      case Ptrac::VALUE:
        nps.m_tally = tmp;
        break;
    }
  }

  if( m_format == Ptrac::BIN_PTRAC ) {
    ReadValue(size2);

    if( size1 != size2 ) {
      throw McnpToolsException( "Failed to read binary PTRAC" );
    }
  }

  hist.m_nps = nps;
  // read the events
  while( (int) next_event_type != Ptrac::LST ) {
    int bnk_type = std::abs(static_cast<int>(next_event_type)) % 1000;
    next_event_type = std::abs(static_cast<int>(next_event_type)) - bnk_type;

    const std::vector<int>& all_data_types = m_event_layouts.at(next_event_type);
    PtracEvent event;
    event.m_type = next_event_type;
    event.m_bnktype = bnk_type;
  
  
    if( m_format == Ptrac::BIN_PTRAC){
      ReadValue(size1);
      binary_buffer.resize(all_data_types.size());
      m_handle.read(reinterpret_cast<char*>(binary_buffer.data()), all_data_types.size() * sizeof(double));
      for (size_t i = 0; i < all_data_types.size(); ++i) {
          double val = binary_buffer[i];
          if (all_data_types[i] == Ptrac::NEXT_EVENT_TYPE) {
              next_event_type = val;
          } else {
              event.m_data.emplace(all_data_types[i], val);
          }
      }
    } 
    else{
        for(unsigned int i=0; i<all_data_types.size(); i++) {
          double tmp;
    
          ReadValue(tmp);
    
          switch( all_data_types[i] ) {
            case Ptrac::NEXT_EVENT_TYPE:
              next_event_type = tmp;
              break;
            case Ptrac::NODE:
            case Ptrac::NSR:
            case Ptrac::ZAID:
            case Ptrac::RXN:
            case Ptrac::SURFACE:
            case Ptrac::ANGLE:
            case Ptrac::TERMINATION_TYPE:
            case Ptrac::BRANCH:
            case Ptrac::PARTICLE:
            case Ptrac::CELL:
            case Ptrac::MATERIAL:
            case Ptrac::COLLISION_NUMBER:
            case Ptrac::X:
            case Ptrac::Y:
            case Ptrac::Z:
            case Ptrac::U:
            case Ptrac::V:
            case Ptrac::W:
            case Ptrac::ENERGY:
            case Ptrac::WEIGHT:
            case Ptrac::TIME:
              event.m_data.emplace(all_data_types[i], tmp);
              break;
          }
        }
    }
    if( m_format == Ptrac::BIN_PTRAC ) {
      ReadValue(size2);

      if( size1 != size2 ) {
        throw McnpToolsException( "Failed to read binary PTRAC" );
      }
    }
    hist.m_events.push_back(std::move(event));
  }

  // make sure to read until end of line if ASCII
  if( m_format == Ptrac::ASC_PTRAC ) {
    std::string chkstr;
    cjsoft::stringops::getline(m_handle, chkstr);
  }

  return hist;
}

std::vector<PtracHistory> Ptrac::ReadHistoriesLegacy(const unsigned int num) {
  std::vector<PtracHistory> retval;
  retval.reserve(num);
  for(unsigned int i=0; i<num; i++) {
    m_handle.peek();
    if( !m_handle.eof() ) {
      retval.emplace_back( ReadHistory() );
    }
    else {
      break;
    }
  }

  return retval;
}


std::vector<PtracHistory> Ptrac::ReadHistories(const unsigned int num) {
  return (m_format == Ptrac::HDF5_PTRAC) ? m_hdf5_parser->ReadHistories( num ) :
                                           ReadHistoriesLegacy( num );
}

} // end namespace mcnptools
