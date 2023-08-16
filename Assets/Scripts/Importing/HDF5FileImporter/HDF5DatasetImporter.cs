using AS.HDFql;
using UnityEngine;
using System;
using System.IO;
using System.Threading.Tasks;

// This is mostly based off of the RawDatasetImporter from the original project!
// Just needed some HDFql code and some array flattening to get there

namespace UnityVolumeRendering
{

    public class HDF5DatasetImporter
    {
        string filePath;
        string dataset;
        private int dimX;
        private int dimY;
        private int dimZ;
        private DataContentFormat contentFormat;
        public HDF5DatasetImporter(string filePath, string dataset, int dimX, int dimY, int dimZ, DataContentFormat contentFormat)
        {
            this.filePath = filePath;
            this.dataset = dataset;
            this.dimX = dimX;
            this.dimY = dimY;
            this.dimZ = dimZ;
            this.contentFormat = contentFormat;
        }

        public VolumeDataset Import()
        {
            // Check that the file exists
            if (!File.Exists(filePath))
            {
                Debug.LogError("The file does not exist: " + filePath);
                return null;
            }

            VolumeDataset volumeDataset = ScriptableObject.CreateInstance<VolumeDataset>();
            ImportInternal(volumeDataset);

            return volumeDataset;
        }
        public async Task<VolumeDataset> ImportAsync()
        {
            // Check that the file exists
            if (!File.Exists(filePath))
            {
                Debug.LogError("The file does not exist: " + filePath);
                return null;
            }

            VolumeDataset volumeDataset = ScriptableObject.CreateInstance<VolumeDataset>();

            await Task.Run(() => ImportInternal(volumeDataset));

            return volumeDataset;
        }
        private void ImportInternal(VolumeDataset volumeDataset)
        {
            volumeDataset.datasetName = Path.GetFileName(filePath);
            volumeDataset.filePath = filePath;
            volumeDataset.dimX = dimX;
            volumeDataset.dimY = dimY;
            volumeDataset.dimZ = dimZ;

            int uDimension = dimX * dimY * dimZ;
            volumeDataset.data = new float[uDimension];

            switch (contentFormat)
            {
                case DataContentFormat.Int8: 
                    {
                        sbyte[,,] temp = pullDataFromFile<sbyte>();
                        for (int i = 0; i < dimX; i++) {
                            for (int j = 0; j < dimY; j++) {
                                for (int k = 0; k < dimZ; k++) {
                                    volumeDataset.data[i+j*dimX+k*dimX*dimY] = (float)temp[i,j,k];
                                }
                            }
                        }                        
                        break;
                    }
                case DataContentFormat.Int16:
                    {
                        short[,,] temp = pullDataFromFile<short>();
                        for (int i = 0; i < dimX; i++) {
                            for (int j = 0; j < dimY; j++) {
                                for (int k = 0; k < dimZ; k++) {
                                    volumeDataset.data[i+j*dimX+k*dimX*dimY] = (float)temp[i,j,k];
                                }
                            }  
                        }                      
                        break;
                    }
                case DataContentFormat.Int32:
                    {
                        int[,,] temp = pullDataFromFile<int>();
                        for (int i = 0; i < dimX; i++) {
                            for (int j = 0; j < dimY; j++) {
                                for (int k = 0; k < dimZ; k++) {
                                    volumeDataset.data[i+j*dimX+k*dimX*dimY] = (float)temp[i,j,k];
                                }
                            }      
                        }                  
                        break;
                    }
                case DataContentFormat.Uint8:
                    {
                        byte[,,] temp = pullDataFromFile<byte>();
                        for (int i = 0; i < dimX; i++) {
                            for (int j = 0; j < dimY; j++) {
                                for (int k = 0; k < dimZ; k++) {
                                    volumeDataset.data[i+j*dimX+k*dimX*dimY] = (float)temp[i,j,k];
                                }
                            }     
                        }                   
                        break;
                    }
                case DataContentFormat.Uint16:
                    {
                        ushort[,,] temp = pullDataFromFile<ushort>();
                        for (int i = 0; i < dimX; i++) {
                            for (int j = 0; j < dimY; j++) {
                                for (int k = 0; k < dimZ; k++) {
                                    volumeDataset.data[i+j*dimX+k*dimX*dimY] = (float)temp[i,j,k];
                                }
                            }     
                        }                   
                        break;
                    }
                case DataContentFormat.Uint32:
                    {
                        uint[,,] temp = pullDataFromFile<uint>();
                        for (int i = 0; i < dimX; i++) {
                            for (int j = 0; j < dimY; j++) {
                                for (int k = 0; k < dimZ; k++) {
                                    volumeDataset.data[i+j*dimX+k*dimX*dimY] = (float)temp[i,j,k];
                                }
                            }
                        }
                        break;
                    }
                default:
                    throw new NotImplementedException("Unimplemented data content format");
            }

            Debug.Log("Loaded dataset in range: " + volumeDataset.GetMinDataValue() + "  -  " + volumeDataset.GetMaxDataValue());

            volumeDataset.FixDimensions();
            volumeDataset.rotation = Quaternion.Euler(90.0f, 0.0f, 0.0f);
        }

        private T[,,] pullDataFromFile<T>() {
            T[,,] temp = new T[dimX, dimY, dimZ];

            // Need to swap to double slashes and include quotes in the filepath so C#/HDFql can read it properly
            HDFql.Execute("USE FILE " + "\"" + filePath.Replace("/","\\") + "\"");

            // Set a variable register (so HDFql knows where temp is in memory) then put the data there
            HDFql.Execute("SELECT FROM " + dataset + " INTO MEMORY " + HDFql.VariableRegister(temp));

            return temp;
        }

        private int GetSampleFormatSize(DataContentFormat format)
        {
            switch (format)
            {
                case DataContentFormat.Int8:
                    return 1;
                case DataContentFormat.Uint8:
                    return 1;
                case DataContentFormat.Int16:
                    return 2;
                case DataContentFormat.Uint16:
                    return 2;
                case DataContentFormat.Int32:
                    return 4;
                case DataContentFormat.Uint32:
                    return 4;
            }
            throw new NotImplementedException();
        }
    }
}
