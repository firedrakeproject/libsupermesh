/*  Copyright (C) 2006 Imperial College London and others.
    
    Please see the AUTHORS file in the main source directory for a full list
    of copyright holders.

    Prof. C Pain
    Applied Modelling and Computation Group
    Department of Earth Science and Engineering
    Imperial College London

    amcgsoftware@imperial.ac.uk
    
    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation,
    version 2.1 of the License.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
    USA
*/

#include "Halos_IO.h"

using namespace std;

using namespace LibSupermesh;

HaloReadError LibSupermesh::ReadHalos(const string& filename, int& process, int& nprocs, map<int, int>& npnodes, map<int, vector<vector<int> > >& send, map<int, vector<vector<int> > >& recv){ 
  // Read the halo file
  TiXmlDocument doc(filename);
  if(!doc.LoadFile()){
    doc.ErrorDesc();
    return HALO_READ_FILE_NOT_FOUND;
  }
  
  const char* charBuffer;
  ostringstream buffer;
   
  // Extract the XML header
  TiXmlNode* header = doc.FirstChild();
  while(header != NULL and header->Type() != TiXmlNode::DECLARATION){
    header = header->NextSibling();
  }

  // Extract the root node
  TiXmlNode* rootNode = header->NextSiblingElement();
  if(rootNode == NULL){
    return HALO_READ_FILE_INVALID;
  }
  TiXmlElement* rootEle = rootNode->ToElement();
  
  // Extract process
  charBuffer = rootEle->Attribute("process");
  if(charBuffer == NULL){
    return HALO_READ_FILE_INVALID;
  }
  process = atoi(charBuffer);
  if(process < 0){
    return HALO_READ_FILE_INVALID;
  }
  
  // Extract nprocs
  charBuffer = rootEle->Attribute("nprocs");
  if(charBuffer == NULL){
    return HALO_READ_FILE_INVALID;
  }
  nprocs = atoi(charBuffer);
  if(process >= nprocs){
    return HALO_READ_FILE_INVALID;
  }
  
  // Extract halo data for each process for each level
  npnodes.clear();
  send.clear();
  recv.clear();
  // Find the next halo element
  for(TiXmlNode* haloNode = rootEle->FirstChildElement("halo");haloNode != NULL;haloNode = haloNode->NextSiblingElement("halo")){
    if(haloNode == NULL){
      break;
    }
    TiXmlElement* haloEle = haloNode->ToElement();
    
    // Extract the level
    charBuffer = haloEle->Attribute("level");
    if(charBuffer == NULL){
      // Backwards compatibility
      charBuffer = haloEle->Attribute("tag");
      if(charBuffer == NULL){
        return HALO_READ_FILE_INVALID;
      }
    }
    int level = atoi(charBuffer);
    send[level] = vector<vector<int> >(nprocs);
    recv[level] = vector<vector<int> >(nprocs);

    // Extract n_private_nodes
    charBuffer = haloEle->Attribute("n_private_nodes");
    if(charBuffer == NULL){
      return HALO_READ_FILE_INVALID;
    }
    npnodes[level] = atoi(charBuffer);
    if(npnodes[level] < 0){
      return HALO_READ_FILE_INVALID;
    }
    
    // Find the next halo_data element
    for(TiXmlNode* dataNode = haloEle->FirstChildElement("halo_data");dataNode != NULL;dataNode = dataNode->NextSiblingElement("halo_data")){
      if(dataNode == NULL){
        break;
      }
      TiXmlElement* dataEle = dataNode->ToElement();
    
      // Extract the process
      charBuffer = dataEle->Attribute("process");
      if(charBuffer == NULL){
        return HALO_READ_FILE_INVALID;
      }
      int proc = atoi(charBuffer);
      if(proc < 0 or proc >= nprocs){
        return HALO_READ_FILE_INVALID;
      }
      
      // Check that data for this level and process has not already been extracted
      if(send[level][proc].size() > 0 or recv[level][proc].size() > 0){
        return HALO_READ_FILE_INVALID;
      }
      
      // Permit empty send and receive data elements
      send[level][proc] = vector<int>();
      recv[level][proc] = vector<int>();
      
      // Extract the send data
      TiXmlNode* sendDataNode = dataEle->FirstChildElement("send");
      if(sendDataNode != NULL){
        TiXmlNode* sendDataTextNode = sendDataNode->FirstChild();
        while(sendDataTextNode != NULL and sendDataTextNode->Type() != TiXmlNode::TEXT){
          sendDataTextNode = sendDataTextNode->NextSibling();
        }
        if(sendDataTextNode != NULL){
          vector<string> tokens;
          Tokenize(sendDataTextNode->ValueStr(), tokens, " ");
          for(size_t i = 0;i < tokens.size();i++){
            send[level][proc].push_back(atoi(tokens[i].c_str()));
          }
        }
      }
      
      // Extract the receive data
      TiXmlNode* recvDataNode = dataEle->FirstChildElement("receive");
      if(recvDataNode != NULL){
      TiXmlNode* recvDataTextNode = recvDataNode->FirstChild();
        while(recvDataTextNode != NULL and recvDataTextNode->Type() != TiXmlNode::TEXT){
          recvDataTextNode = recvDataTextNode->NextSibling();
        }
        if(recvDataTextNode != NULL){
          vector<string> tokens;
          Tokenize(recvDataTextNode->ValueStr(), tokens, " ");
          for(size_t i = 0;i < tokens.size();i++){
            recv[level][proc].push_back(atoi(tokens[i].c_str()));
          }
        }
      }
    }
  }

  return HALO_READ_SUCCESS;
}

namespace LibSuperMesh{
  HaloData* readHaloData = NULL;
}
  
extern "C"{
  void cLibSuperMesh_halo_reader_reset(){
    if(LibSuperMesh::readHaloData){
      delete LibSuperMesh::readHaloData;
      LibSuperMesh::readHaloData = NULL;
    }
  
    return;
  }
  
  int cLibSuperMesh_halo_reader_set_input(char* filename, int* filename_len, int* process, int* nprocs){    
    if(LibSuperMesh::readHaloData){
      delete LibSuperMesh::readHaloData;
      LibSuperMesh::readHaloData = NULL;
    }
    LibSuperMesh::readHaloData = new HaloData();    
  
    ostringstream buffer;
    buffer << string(filename, *filename_len) << "_" << *process << ".halo";
    HaloReadError ret = ReadHalos(buffer.str(), 
      LibSuperMesh::readHaloData->process, LibSuperMesh::readHaloData->nprocs,
      LibSuperMesh::readHaloData->npnodes, LibSuperMesh::readHaloData->send, LibSuperMesh::readHaloData->recv);
  
    int errorCount = 0;
    if(ret == HALO_READ_FILE_NOT_FOUND){
      if(*process == 0){
        cerr << "Error reading halo file " << buffer.str() << "\n"
             << "Zero process file not found" << endl;
        errorCount++;
      }else{
        LibSuperMesh::readHaloData->process = *process;
        LibSuperMesh::readHaloData->nprocs = *nprocs;
        LibSuperMesh::readHaloData->npnodes.clear();
        LibSuperMesh::readHaloData->send.clear();
        LibSuperMesh::readHaloData->recv.clear();
      }
    }else if(ret != HALO_READ_SUCCESS){
      cerr << "Error reading halo file " << buffer.str() << "\n";
      switch(ret){
        case(HALO_READ_FILE_INVALID):
          cerr << "Invalid .halo file" << endl;
          break;
        // HALO_READ_FILE_NOT_FOUND case handled above
        default:
          cerr << "Unknown error" << endl;
          break;
      }
      errorCount++;
    }else if(LibSuperMesh::readHaloData->process != *process){
      cerr << "Error reading halo file " << buffer.str() << "\n"
           << "Unexpected process number in .halo file" << endl;
      errorCount++;
    }else if(LibSuperMesh::readHaloData->nprocs > *nprocs){
      cerr << "Error reading halo file " << buffer.str() << "\n"
           << "Number of processes in .halo file: " << LibSuperMesh::readHaloData->nprocs << "\n"
           << "Number of running processes: " << *nprocs << "\n"
           << "Number of processes in .halo file exceeds number of running processes" << endl;
      errorCount++;
    }
    
    return errorCount;
  }
  
  void cLibSuperMesh_halo_reader_query_output(int* level, int* nprocs, int* nsends, int* nreceives){
    assert(LibSuperMesh::readHaloData);
    assert(*nprocs >= LibSuperMesh::readHaloData->nprocs);
    
    if(LibSuperMesh::readHaloData->send.count(*level) == 0){
      assert(LibSuperMesh::readHaloData->recv.count(*level) == 0);
      
      for(int i = 0;i < *nprocs;i++){
        nsends[i] = 0;
        nreceives[i] = 0;
      }
    }else{
      assert(LibSuperMesh::readHaloData->recv.count(*level) > 0);
      
      for(int i = 0;i < LibSuperMesh::readHaloData->nprocs;i++){
        nsends[i] = LibSuperMesh::readHaloData->send[*level][i].size();
        nreceives[i] = LibSuperMesh::readHaloData->recv[*level][i].size();
      }
      for(int i = LibSuperMesh::readHaloData->nprocs;i < *nprocs;i++){
        nsends[i] = 0;
        nreceives[i] = 0;
      }
    }
    
    return;
  }
  
  void cLibSuperMesh_halo_reader_get_output(int* level, int* nprocs, int* nsends, int* nreceives,
    int* npnodes, int* send, int* recv){
    
#ifdef DDEBUG
    assert(LibSuperMesh::readHaloData);
    assert(*nprocs >= LibSuperMesh::readHaloData->nprocs);
    int* lnsends = (int*)malloc(*nprocs * sizeof(int));
    assert(lnsends);
    int* lnreceives = (int*)malloc(*nprocs * sizeof(int));
    assert(lnreceives);
    cLibSuperMesh_halo_reader_query_output(level, nprocs, lnsends, lnreceives);
    for(int i = 0;i < *nprocs;i++){
      assert(nsends[i] == lnsends[i]);
      assert(nreceives[i] == lnreceives[i]);
    }
    free(lnsends);
    free(lnreceives);
#endif

    if(LibSuperMesh::readHaloData->send.count(*level) == 0){
#ifdef DDEBUG
      assert(LibSuperMesh::readHaloData->recv.count(*level) == 0);
      for(int i = 0;i < *nprocs;i++){
        assert(nsends[i] == 0);
        assert(nreceives[i] == 0);
      }
#endif
    }else{
      assert(LibSuperMesh::readHaloData->recv.count(*level) > 0);
      
      int sendIndex = 0, recvIndex = 0;;
      for(int i = 0;i < LibSuperMesh::readHaloData->nprocs;i++){
        for(int j = 0;j < nsends[i];j++){
          send[sendIndex] = LibSuperMesh::readHaloData->send[*level][i][j];
          sendIndex++;
        }
        for(int j = 0;j < nreceives[i];j++){
          recv[recvIndex] = LibSuperMesh::readHaloData->recv[*level][i][j];
          recvIndex++;
        }
      }
    }
    
    *npnodes = LibSuperMesh::readHaloData->npnodes[*level];
    
    return;
    }
}
