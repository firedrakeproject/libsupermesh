/*
 * The following code is derived from femtools/Halos_IO.cpp and
 * femtools/Tokenize.cpp in Fluidity git
 * revision 4e6c1d2b022df3a519cdec120fad28e60d1b08d9 (dated 2015-02-25)
 */
 
// Fluidity copyright information (note that AUTHORS mentioned in the following
// has been renamed to fluidity_AUTHORS):

/*
    Copyright (C) 2006 Imperial College London and others.
    
    Please see the AUTHORS file in the main source directory for a full list
    of copyright holders.

    Prof. C Pain
    Applied Modelling and Computation Group
    Department of Earth Science and Engineering
    Imperial College London

    amcgsoftware@imperial.ac.uk
*/

#include "Read_Halos_C++.h"

using namespace std;
using namespace libsupermesh;

void libsupermesh::ReadHalos(const string &filename, int &process, int &nprocs,
  map<int, int> &npnodes, map<int, vector<vector<int> > > &send,
  map<int, vector<vector<int> > > &recv) { 
  npnodes.clear();
  send.clear();
  recv.clear();
  
  // Read the halo file
  TiXmlDocument doc(filename.c_str());
  if(!doc.LoadFile()) {
    cerr << doc.ErrorDesc() << endl;
    cerr << "Failed to read Halo file '" << filename << "'" << endl;
    exit(1);
  }
  
  const char *buffer;
   
  // Extract the XML header
  TiXmlNode *header = doc.FirstChild();
  while(header and header->Type() != TiXmlNode::TINYXML_DECLARATION) {
    header = header->NextSibling();
  }
  if(!header) {
    cerr << "Invalid halo file '" << filename << "': Missing XML declaration" << endl;
    exit(1);
  }

  // Extract the root node
  TiXmlNode *rootNode = header->NextSiblingElement();
  if(!rootNode) {
    cerr << "Invalid halo file '" << filename << "': Missing root element" << endl;
    exit(1);
  }
  TiXmlElement *rootEle = rootNode->ToElement();
  
  // Extract process
  buffer = rootEle->Attribute("process");
  if(!buffer) {
    cerr << "Invalid halo file '" << filename << "': Missing process attribute" << endl;
    exit(1);
  }
  process = atoi(buffer);
  if(process < 0) {
    cerr << "Invalid halo file '" << filename << "': Invalid process attribute" << endl;
    exit(1);
  }
  
  // Extract nprocs
  buffer = rootEle->Attribute("nprocs");
  if(!buffer) {
    cerr << "Invalid halo file '" << filename << "': Missing nprocs attribute" << endl;
    exit(1);
  }
  nprocs = atoi(buffer);
  if(nprocs < 1) {
    cerr << "Invalid halo file '" << filename << "': Invalid nprocs attribute" << endl;
    exit(1);
  } else if(process >= nprocs) {
    cerr << "Invalid halo file '" << filename << "': Invalid process / nprocs attributes" << endl;
    exit(1);
  }
  
  // Extract halo data for each process for each level
  // Find the next halo element
  for(TiXmlNode *haloNode = rootEle->FirstChildElement("halo");haloNode;haloNode = haloNode->NextSiblingElement("halo")) {
    TiXmlElement *haloEle = haloNode->ToElement();
    
    // Extract the level
    buffer = haloEle->Attribute("level");
    if(!buffer) {
      cerr << "Invalid halo file '" << filename << "': halo_data element missing level attribute" << endl;
      exit(1);
    }
    // Check that data for this level has not already been extracted
    int level = atoi(buffer);
    if(send.count(level) > 0 or recv.count(level) > 0) {
      cerr << "Invalid halo file '" << filename << "': Multiple halos defined for a single level" << endl;
      exit(1);
    }
    send[level] = vector<vector<int> >(nprocs);
    recv[level] = vector<vector<int> >(nprocs);

    // Extract n_private_nodes
    buffer = haloEle->Attribute("n_private_nodes");
    if(!buffer) {      
      cerr << "Invalid halo file '" << filename << "': halo_data element missing n_private_nodes attribute" << endl;
      exit(1);
    }
    npnodes[level] = atoi(buffer);
    
    // Find the next halo_data element
    for(TiXmlNode *dataNode = haloEle->FirstChildElement("halo_data");dataNode;dataNode = dataNode->NextSiblingElement("halo_data")) {
      TiXmlElement *dataEle = dataNode->ToElement();
    
      // Extract the process
      buffer = dataEle->Attribute("process");
      if(!buffer) {
        cerr << "Invalid halo file '" << filename << "': halo_data element missing process attribute" << endl;
        exit(1);
      }
      int proc = atoi(buffer);
      if(proc < 0 or proc >= nprocs) {
        cerr << "Invalid halo file '" << filename << "': Invalid halo_data element process / nprocs attributes" << endl;
        exit(1);
      }
      
      // Check that data for this level and process has not already been extracted
      if(send[level][proc].size() > 0 or recv[level][proc].size() > 0) {        
        cerr << "Invalid halo file '" << filename << "': Multiple halos defined for a single level and process" << endl;
        exit(1);
      }
      
      // Extract the send data
      TiXmlNode *sendDataNode = dataEle->FirstChildElement("send");
      if(sendDataNode) {
        TiXmlNode *sendDataTextNode = sendDataNode->FirstChild();
        while(sendDataTextNode and sendDataTextNode->Type() != TiXmlNode::TINYXML_TEXT) {
          sendDataTextNode = sendDataTextNode->NextSibling();
        }
        if(sendDataTextNode) {
          vector<string> tokens;
          Tokenize(string(sendDataTextNode->Value()), tokens, " ");
          for(vector<string>::size_type i = 0;i < tokens.size();i++) {
            send[level][proc].push_back(atoi(tokens[i].c_str()));
          }
        }
        
        if(sendDataNode->NextSiblingElement("data")) {        
          cerr << "Invalid halo file '" << filename << "': Multiple sets of sends defined for a single level and process" << endl;
          exit(1);
        }
      }
      
      // Extract the receive data
      TiXmlNode *recvDataNode = dataEle->FirstChildElement("receive");
      if(recvDataNode) {
        TiXmlNode *recvDataTextNode = recvDataNode->FirstChild();
        while(recvDataTextNode and recvDataTextNode->Type() != TiXmlNode::TINYXML_TEXT) {
          recvDataTextNode = recvDataTextNode->NextSibling();
        }
        if(recvDataTextNode) {
          vector<string> tokens;
          Tokenize(string(recvDataTextNode->Value()), tokens, " ");
          for(vector<string>::size_type i = 0;i < tokens.size();i++) {
            recv[level][proc].push_back(atoi(tokens[i].c_str()));
          }
        }
        
        if(recvDataNode->NextSiblingElement("data")) {        
          cerr << "Invalid halo file '" << filename << "': Multiple sets of receives defined for a single level and process" << endl;
          exit(1);
        }
      }
    }
  }

  return;
}

void libsupermesh::Tokenize(const string &str, vector<string> &tokens, const string &delimiters) {
  tokens.clear();
  
  // Skip delimiter at beginning
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  
  // Find first delimiter
  string::size_type pos = str.find_first_of(delimiters, lastPos);
  
  while(lastPos != string::npos) {
    // Found a token, add it to the vector
    tokens.push_back(str.substr(lastPos, pos - lastPos));

    // Skip delimiters.  Note the "not_of".
    lastPos = str.find_first_not_of(delimiters, pos);

    // Find next delimiter
    pos = str.find_first_of(delimiters, lastPos);
  }
  
  return;
}
  
extern "C" {  
  void libsupermesh_read_halo(void **data, const char *basename,
    const int *basename_len, const int *process, const int *nprocs) {      
    (*data) = (void*)(new HaloData());
    if(!((HaloData*)(*data))) {
      cerr << "new failure" << endl;
      exit(1);
    }
    
    ostringstream filename;
    filename << string(basename, *basename_len) << "_" << *process << ".halo";
  
    ReadHalos(filename.str(),
      ((HaloData*)(*data))->process, ((HaloData*)(*data))->nprocs,
      ((HaloData*)(*data))->npnodes, ((HaloData*)(*data))->send, ((HaloData*)(*data))->recv);  
    if(((HaloData*)(*data))->process != *process) {    
      cerr << "Failed to read Halo file '" << filename << "': Unexpected process number" << endl;
      exit(1);
    } else if(((HaloData*)(*data))->nprocs != *nprocs) {    
      cerr << "Failed to read Halo file '" << filename << "': Unexpected number of processes" << endl;
      exit(1);
    }
    
    return;
  }
  
  void libsupermesh_halo_sizes(const void **data, const int *level,
    const int *nprocs, int *nsends, int *nreceives) {
    assert(*nprocs == ((HaloData*)(*data))->nprocs);
    
    if(((HaloData*)(*data))->send.count(*level) == 0) {      
      for(int i = 0;i < *nprocs;i++) {
        nsends[i] = 0;
      }
    } else {      
      for(int i = 0;i < *nprocs;i++) {
        nsends[i] = ((HaloData*)(*data))->send[*level][i].size();
      }
    }
    if(((HaloData*)(*data))->recv.count(*level) == 0) {      
      for(int i = 0;i < *nprocs;i++) {
        nreceives[i] = 0;
      }
    } else {      
      for(int i = 0;i < *nprocs;i++) {
        nreceives[i] = ((HaloData*)(*data))->recv[*level][i].size();
      }
    }
    
    return;
  }
  
  void libsupermesh_halo_data(const void **data, const int *level,
    const int *nprocs, const int *nsends, const int *nreceives, int *npnodes,
    int *send, int *recv) {
#ifndef NDEBUG
    assert(*nprocs == ((HaloData*)(*data))->nprocs);
    for(int i = 0;i < *nprocs;i++) {
      assert(nsends[i] == ((HaloData*)(*data))->send[*level][i].size());
      assert(nreceives[i] == ((HaloData*)(*data))->recv[*level][i].size());
    }
#endif
    
    *npnodes = ((HaloData*)(*data))->npnodes[*level];
    if(((HaloData*)(*data))->send.count(*level) > 0) {
      int sendIndex = 0;
      for(int i = 0;i < *nprocs;i++) {
        for(int j = 0;j < nsends[i];j++) {
          send[sendIndex] = ((HaloData*)(*data))->send[*level][i][j];
          sendIndex++;
        }
      }
    }
    if(((HaloData*)(*data))->recv.count(*level) > 0) {
      int recvIndex = 0;
      for(int i = 0;i < *nprocs;i++) {
        for(int j = 0;j < nreceives[i];j++) {
          recv[recvIndex] = ((HaloData*)(*data))->recv[*level][i][j];
          recvIndex++;
        }
      }
    }
    
    return;
    }
    
  void libsupermesh_deallocate_halo(void **data) {
    delete ((HaloData*)(*data));
  
    return;
  }
}
