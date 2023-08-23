using AS.HDFql;
using UnityEngine;
using System;
using System.IO;
using System.Threading.Tasks;

namespace UnityVolumeRendering
{
    [Serializable]
    public enum CoordinateSystem {
        Cartesian,
        Spherical,
    }

    [Serializable]
    public enum AngleUnits {
        Radians,
        Degrees
    }

    public class HDF5DatasetImporter
    {
        string filePath;
        string dataset;
        private int[] dataSize;
        private DataContentFormat contentFormat;
        private CoordinateSystem coordinateSystem;
        private AngleUnits angleUnits;
        private float rMin;
        private float rMax;
        private float thetaMin;
        private float thetaMax;
        private float phiMin;
        private float phiMax;
        private int[] gridSize;
        private bool filterData;
        private float filterValue;

        string rData;
        string thetaData;
        string phiData;

        public HDF5DatasetImporter() 
        {

        }

        public HDF5DatasetImporter(string filePath, string dataset, int[] dataSize, DataContentFormat contentFormat,
                                    CoordinateSystem coordinateSystem)
        {
            this.filePath = filePath;
            this.dataset = dataset;
            this.dataSize = dataSize;
            this.contentFormat = contentFormat;
            this.coordinateSystem = coordinateSystem;            
        }

        public HDF5DatasetImporter(string filePath, string dataset, int[] dataSize, DataContentFormat contentFormat,
                                   CoordinateSystem coordinateSystem, AngleUnits angleUnits, float rMin, float rMax, float thetaMin,
                                   float thetaMax, float phiMin, float phiMax, int[] gridSize, bool filterData, float filterValue)
        {
            this.filePath = filePath;
            this.dataset = dataset;
            this.dataSize = dataSize;
            this.contentFormat = contentFormat;
            this.coordinateSystem = coordinateSystem;
            this.angleUnits = angleUnits;
            this.rMin = rMin;
            this.rMax = rMax;
            this.thetaMin = thetaMin;
            this.thetaMax = thetaMax;
            this.phiMin = phiMin;
            this.phiMax = phiMax;
            this.gridSize = gridSize;
            this.filterData = filterData;
            this.filterValue = filterValue;
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

            if (coordinateSystem == CoordinateSystem.Spherical) {
                volumeDataset.dimX = gridSize[0];
                volumeDataset.dimY = gridSize[1];
                volumeDataset.dimZ = gridSize[2];
            } else {
                volumeDataset.dimX = dataSize[0];
                volumeDataset.dimY = dataSize[1];
                volumeDataset.dimZ = dataSize[2];
            }

            switch (contentFormat)
            {
                case DataContentFormat.Int8: 
                    {
                        // get the data out of the HDF5 file
                        sbyte[,,] temp = pullDataFromFile<sbyte>();

                        // convert it to floats
                        float[,,] data = new float[dataSize[0], dataSize[1], dataSize[2]];

                        for (int i = 0; i < dataSize[0]; i++) {
                            for (int j = 0; j < dataSize[1]; j++) {
                                for (int k = 0; k < dataSize[2]; k++) {
                                    data[i,j,k] = (float)temp[i,j,k];
                                }
                            }
                        }

                        if (coordinateSystem == CoordinateSystem.Spherical) {
                            data = sphericalToCartesianGrid(data);

                        }
                        volumeDataset.data = flattenData(data, volumeDataset.dimX, volumeDataset.dimY, volumeDataset.dimZ);                       
                        break;
                    }
                case DataContentFormat.Int16:
                    {
                        short[,,] temp = pullDataFromFile<short>();

                        // convert it to floats
                        float[,,] data = new float[dataSize[0], dataSize[1], dataSize[2]];

                        for (int i = 0; i < dataSize[0]; i++) {
                            for (int j = 0; j < dataSize[1]; j++) {
                                for (int k = 0; k < dataSize[2]; k++) {
                                    data[i,j,k] = (float)temp[i,j,k];
                                }
                            }
                        }

                        if (coordinateSystem == CoordinateSystem.Spherical) {
                            data = sphericalToCartesianGrid(data);

                        }
                        volumeDataset.data = flattenData(data, volumeDataset.dimX, volumeDataset.dimY, volumeDataset.dimZ);                     
                        break;
                    }
                case DataContentFormat.Int32:
                    {
                        int[,,] temp = pullDataFromFile<int>();

                        // convert it to floats
                        float[,,] data = new float[dataSize[0], dataSize[1], dataSize[2]];

                        for (int i = 0; i < dataSize[0]; i++) {
                            for (int j = 0; j < dataSize[1]; j++) {
                                for (int k = 0; k < dataSize[2]; k++) {
                                    data[i,j,k] = (float)temp[i,j,k];
                                }
                            }
                        }

                        if (coordinateSystem == CoordinateSystem.Spherical) {
                            data = sphericalToCartesianGrid(data);

                        }
                        volumeDataset.data = flattenData(data, volumeDataset.dimX, volumeDataset.dimY, volumeDataset.dimZ);               
                        break;
                    }
                case DataContentFormat.Uint8:
                    {
                        byte[,,] temp = pullDataFromFile<byte>();

                        // convert it to floats
                        float[,,] data = new float[dataSize[0], dataSize[1], dataSize[2]];

                        for (int i = 0; i < dataSize[0]; i++) {
                            for (int j = 0; j < dataSize[1]; j++) {
                                for (int k = 0; k < dataSize[2]; k++) {
                                    data[i,j,k] = (float)temp[i,j,k];
                                }
                            }
                        }

                        if (coordinateSystem == CoordinateSystem.Spherical) {
                            data = sphericalToCartesianGrid(data);

                        }
                        volumeDataset.data = flattenData(data, volumeDataset.dimX, volumeDataset.dimY, volumeDataset.dimZ);                
                        break;
                    }
                case DataContentFormat.Uint16:
                    {
                        ushort[,,] temp = pullDataFromFile<ushort>();

                        // convert it to floats
                        float[,,] data = new float[dataSize[0], dataSize[1], dataSize[2]];

                        for (int i = 0; i < dataSize[0]; i++) {
                            for (int j = 0; j < dataSize[1]; j++) {
                                for (int k = 0; k < dataSize[2]; k++) {
                                    data[i,j,k] = (float)temp[i,j,k];
                                }
                            }
                        }

                        if (coordinateSystem == CoordinateSystem.Spherical) {
                            data = sphericalToCartesianGrid(data);

                        }
                        volumeDataset.data = flattenData(data, volumeDataset.dimX, volumeDataset.dimY, volumeDataset.dimZ);            
                        break;
                    }
                case DataContentFormat.Uint32:
                    {
                        uint[,,] temp = pullDataFromFile<uint>();

                        // convert it to floats
                        float[,,] data = new float[dataSize[0], dataSize[1], dataSize[2]];

                        for (int i = 0; i < dataSize[0]; i++) {
                            for (int j = 0; j < dataSize[1]; j++) {
                                for (int k = 0; k < dataSize[2]; k++) {
                                    data[i,j,k] = (float)temp[i,j,k];
                                }
                            }
                        }

                        if (coordinateSystem == CoordinateSystem.Spherical) {
                            data = sphericalToCartesianGrid(data);

                        }
                        volumeDataset.data = flattenData(data, volumeDataset.dimX, volumeDataset.dimY, volumeDataset.dimZ);
                        break;
                    }
                case DataContentFormat.Float32:
                    {
                        // get the data from the HDF5 file
                        float[,,] temp = pullDataFromFile<float>();

                        // convert it to floats
                        float[,,] data = new float[dataSize[0], dataSize[1], dataSize[2]];

                        for (int i = 0; i < dataSize[0]; i++) {
                            for (int j = 0; j < dataSize[1]; j++) {
                                for (int k = 0; k < dataSize[2]; k++) {
                                    data[i,j,k] = (float)temp[i,j,k];
                                }
                            }
                        }

                        if (coordinateSystem == CoordinateSystem.Spherical) {
                            data = sphericalToCartesianGrid(data);

                        }
                        volumeDataset.data = flattenData(data, volumeDataset.dimX, volumeDataset.dimY, volumeDataset.dimZ);
                        break;
                    }
                case DataContentFormat.Float64:
                    {
                        // get the data from the HDF5 file
                        double[,,] temp = pullDataFromFile<double>();

                        // convert it to floats
                        float[,,] data = new float[dataSize[0], dataSize[1], dataSize[2]];

                        for (int i = 0; i < dataSize[0]; i++) {
                            for (int j = 0; j < dataSize[1]; j++) {
                                for (int k = 0; k < dataSize[2]; k++) {
                                    data[i,j,k] = (float)temp[i,j,k];
                                }
                            }
                        }

                        if (coordinateSystem == CoordinateSystem.Spherical) {
                            data = sphericalToCartesianGrid(data);

                        }
                        volumeDataset.data = flattenData(data, volumeDataset.dimX, volumeDataset.dimY, volumeDataset.dimZ);
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
            T[,,] temp = new T[dataSize[0], dataSize[1], dataSize[2]];

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

        private float[] flattenData(float[,,] data, int firstDim, int secondDim, int thirdDim) {
            float[] flattened = new float[firstDim*secondDim*thirdDim];
            for (int i = 0; i < firstDim; i++) {
                for (int j = 0; j < secondDim; j++) {
                    for (int k  = 0; k < thirdDim; k++) {
                        flattened[i+j*firstDim+k*firstDim*secondDim] = data[i,j,k];
                    }
                }
            }
            return flattened;
        }

        private float[,,] sphericalToCartesianGrid(float[,,] densities) {

            // compute and store r, theta, and phi values, assuming they're equally dispersed across the interval
            float[] rValues = new float[dataSize[0]];
            float[] thetaValues = new float[dataSize[1]];
            float[] phiValues = new float[dataSize[2]];


            for (int i = 0; i < dataSize[0]; i++) {
                rValues[i] = Mathf.Lerp(rMin, rMax, (float)i/(dataSize[0]-1));
            }

            bool degrees = angleUnits == AngleUnits.Degrees;

            for (int i = 0; i < dataSize[1]; i++) {
                if (degrees) {
                    thetaValues[i] = Mathf.Lerp(thetaMin*(Mathf.PI/180), thetaMax*(Mathf.PI/180), (float)i/(dataSize[1]-1));
                } else {
                    thetaValues[i] = Mathf.Lerp(thetaMin, thetaMax, (float)i/(dataSize[1]-1));
                } 
            }

            for (int i = 0; i < dataSize[2]; i++) {
                if (degrees) {  
                    phiValues[i] = Mathf.Lerp(phiMin*(Mathf.PI/180), phiMax*(Mathf.PI/180), (float)i/(dataSize[2]-1));
                } else {
                    phiValues[i] = Mathf.Lerp(phiMin, phiMax, (float)i/(dataSize[2]-1));
                }
            }

            // use the ranges to build the Cartesian point positions
            // we also flatten the density to make it one-to-one with shape of the positions

            float[] x = new float[dataSize[0]*dataSize[1]*dataSize[2]];
            float[] y = new float[dataSize[0]*dataSize[1]*dataSize[2]];
            float[] z = new float[dataSize[0]*dataSize[1]*dataSize[2]];
            float[] densitiesFlattened = new float[dataSize[0]*dataSize[1]*dataSize[2]];

            for (int i = 0; i < dataSize[0]; i++) {
                for (int j = 0; j < dataSize[1]; j++) {
                    for (int k = 0; k < dataSize[2]; k++) {
                        x[i+j*dataSize[0]+k*dataSize[0]*dataSize[1]] = rValues[i]*Mathf.Sin(thetaValues[j])*Mathf.Cos(phiValues[k]);
                        y[i+j*dataSize[0]+k*dataSize[0]*dataSize[1]] = rValues[i]*Mathf.Sin(thetaValues[j])*Mathf.Sin(phiValues[k]);
                        z[i+j*dataSize[0]+k*dataSize[0]*dataSize[1]] = rValues[i]*Mathf.Cos(thetaValues[j]);
                        densitiesFlattened[i+j*dataSize[0]+k*dataSize[0]*dataSize[1]] = densities[i,j,k];
                    }
                }
            }

            // filter the data if desired
            // this improves the density grid resolution in the area you care about if there are sparse sections of the data 
            if (filterData) {
                ulong count = 0;
                for (int i = 0; i < dataSize[0]*dataSize[1]*dataSize[2]; i++) {
                    if (densitiesFlattened[i] >= filterValue) {
                        count++;
                    }
                }
                float[] xFiltered = new float[count];
                float[] yFiltered = new float[count];
                float[] zFiltered = new float[count];
                float[] densitiesFiltered = new float[count];

                Debug.Log(count);
                Debug.Log(dataSize[0]*dataSize[1]*dataSize[2]);

                ulong index = 0;
                for (int i = 0; i < dataSize[0]*dataSize[1]*dataSize[2]; i++) {
                    if (densitiesFlattened[i] >= filterValue) {
                        xFiltered[index] = x[i];
                        yFiltered[index] = y[i];
                        zFiltered[index] = z[i];
                        densitiesFiltered[index] = densitiesFlattened[i];
                        index++;
                    }
                }

                x = xFiltered;
                y = yFiltered;
                z = zFiltered;
                densitiesFlattened = densitiesFiltered;

            }

            Debug.Log(x.GetLength(0));

            // find the right size for the Cartesian grid cells

            float xMin = get1DMin(x);
            float xMax = get1DMax(x);
            float yMin = get1DMin(y);
            float yMax = get1DMax(y);
            float zMin = get1DMin(z);
            float zMax = get1DMax(z);

            // make each cell a little bit larger than it should be (with the 1.001f)
            // this makes sure floating point error in computing xInd,yInd,zInd doesn't make us index out of the array when x ~ xMax
            float cellX = ((xMax - xMin)/gridSize[0])*1.00001f;
            float cellY = ((yMax - yMin)/gridSize[1])*1.00001f;
            float cellZ = ((zMax - zMin)/gridSize[2])*1.00001f;

            Debug.Log(xMax-xMin);
            Debug.Log(yMax-yMin);
            Debug.Log(zMax-zMin);

            float[,,] densitiesCartesian = new float[gridSize[0],gridSize[1],gridSize[2]];
            int[,,] counts = new int[gridSize[0],gridSize[1],gridSize[2]];

            for (int i = 0 ; i < densitiesFlattened.GetLength(0); i++) {
                ushort xInd = (ushort)Mathf.Floor((x[i]-xMin)/cellX);
                ushort yInd = (ushort)Mathf.Floor((y[i]-yMin)/cellY);
                ushort zInd = (ushort)Mathf.Floor((z[i]-zMin)/cellZ);

                densitiesCartesian[xInd,yInd,zInd] = (densitiesCartesian[xInd,yInd,zInd]*counts[xInd,yInd,zInd] +
                                                    densitiesFlattened[i])/(counts[xInd,yInd,zInd]+1);
                counts[xInd,yInd,zInd]++;
            }

            return densitiesCartesian;
        }

        private float get3DMax(float[,,] data) {
            float max = Single.MinValue;
            for (int i = 0; i < data.GetLength(0); i++) {
                for (int j = 0; j < data.GetLength(1); j++) {
                    for (int k = 0; k < data.GetLength(2); k++) {
                        float value = data[i,j,k];
                        if (value > max) {
                            max = value;
                        }
                    }
                }
            }
            return max;
        }

        private float get3DMin(float[,,] data) {
            float min = Single.MaxValue;
            for (int i = 0; i < data.GetLength(0); i++) {
                for (int j = 0; j < data.GetLength(1); j++) {
                    for (int k = 0; k < data.GetLength(2); k++) {
                        float value = data[i,j,k];
                        if (value < min) {
                            min = value;
                        }
                    }
                }
            }
            return min;
        }

        private float get1DMax(float[] data) {
            float max = Single.MinValue;
            for (int i = 0; i < data.GetLength(0); i++) {
                if (data[i] > max) {
                    max = data[i];
                }
            }
            return max;
        }

        private float get1DMin(float[] data) {
            float min = Single.MaxValue;
            for (int i = 0; i < data.GetLength(0); i++) {
                if (data[i] < min) {
                    min = data[i];
                }
            }
            return min;            
        }
    }
}
