using UnityEngine;
using System;
using System.IO;
using System.Threading.Tasks;
using HDF.PInvoke;
using System.Threading;
using System.Runtime.InteropServices;

namespace UnityVolumeRendering
{
    [Serializable]
    public enum CoordinateSystem {
        Cartesian,
        Spherical,
    }

    [Serializable]
    public enum SphericalType {
        Uniform,
        NonUniform
    }

    [Serializable]
    public enum AngleUnits {
        Radians,
        Degrees
    }

    [Serializable]
    public enum SimulationType {
        GridBased,
        ParticleBased
    }

    public class HDF5DatasetImporter
    {
        // Used in all situations
        string filePath;
        string dataset;
        private int[] dataSize = new int[3];
        private int[] gridSize;
        private DataContentFormat contentFormat = DataContentFormat.Uint8;
        private CoordinateSystem coordinateSystem;
        private SimulationType simType;
        private SphericalType sphericalType;
        private AngleUnits angleUnits;
        private bool filterData;
        private float filterValue;

        // For spherical, uniform grid
        private float rMin;
        private float rMax;
        private float thetaMin;
        private float thetaMax;
        private float phiMin;
        private float phiMax;

        // For spherical, non-uniform grid and spherical particle-based
        private string rData;
        private string thetaData;
        private string phiData;

        // For cartesian particle-based
        private string xData;
        private string yData;
        private string zData;

        public HDF5DatasetImporter() {}

        public HDF5DatasetImporter(string filePath, string dataset, //int[] dataSize, DataContentFormat contentFormat,
                                    SimulationType simType, CoordinateSystem coordinateSystem, bool filterData, float filterValue)
        {
            this.filePath = filePath;
            this.dataset = dataset;
            this.simType = simType;
            this.coordinateSystem = coordinateSystem;      
            this.filterData = filterData;
            this.filterValue = filterValue;   
        }

        public HDF5DatasetImporter(string filePath, string dataset, //int[] dataSize, DataContentFormat contentFormat,
                                   SimulationType simType, CoordinateSystem coordinateSystem, SphericalType sphericalType,
                                   AngleUnits angleUnits, float rMin, float rMax, float thetaMin, float thetaMax, 
                                   float phiMin, float phiMax, int[] gridSize, bool filterData, float filterValue)
        {
            this.filePath = filePath;
            this.dataset = dataset;
            this.simType = simType;
            this.coordinateSystem = coordinateSystem;
            this.sphericalType = sphericalType;
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

        public HDF5DatasetImporter(string filePath, string dataset, //int[] dataSize, DataContentFormat contentFormat,
                                   SimulationType simType, CoordinateSystem coordinateSystem, SphericalType sphericalType,
                                   AngleUnits angleUnits, string rData, string thetaData, string phiData, int[] gridSize, 
                                   bool filterData, float filterValue)
        {
            this.filePath = filePath;
            this.dataset = dataset;
            this.simType = simType;
            this.coordinateSystem = coordinateSystem;
            this.angleUnits = angleUnits;
            this.sphericalType = sphericalType;
            this.rData = rData;
            this.thetaData = thetaData;
            this.phiData = phiData;
            this.gridSize = gridSize;
            this.filterData = filterData;
            this.filterValue = filterValue;
        }


        public HDF5DatasetImporter(string filePath, string dataset, //int[] dataSize, DataContentFormat contentFormat,
                                   SimulationType simType, CoordinateSystem coordinateSystem, AngleUnits angleUnits, 
                                   string xData, string yData, string zData, int[] gridSize, bool filterData, float filterValue)
        {
            this.filePath = filePath;
            this.dataset = dataset;
            this.coordinateSystem = coordinateSystem;
            this.angleUnits = angleUnits;
            this.xData = xData;
            this.yData = yData;
            this.zData = zData;
            this.gridSize = gridSize;
            this.filterData = filterData;
            this.filterValue = filterValue;
        }

        public HDF5DatasetImporter(string filePath, string dataset, //int[] dataSize, DataContentFormat contentFormat,
                                   SimulationType simType, CoordinateSystem coordinateSystem, string rData, string thetaData, 
                                   string phiData, AngleUnits angleUnits, int[] gridSize, bool filterData, float filterValue)
        {
            this.filePath = filePath;
            this.dataset = dataset;
            this.coordinateSystem = coordinateSystem;
            this.angleUnits = angleUnits;
            this.rData = rData;
            this.thetaData = thetaData;  
            this.phiData = phiData;
            this.angleUnits = angleUnits;
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

            Debug.Log(dataSize[0]);
            getDataPropertiesFromFile(dataset);
            Debug.Log(dataSize[0]);

            // We can use the raw data size only when we don't have to bin the data,
            // which is only in grid-based cartesian simulations
            // Otherwise the final data takes the form of the desired grid
            if (simType == SimulationType.GridBased && coordinateSystem == CoordinateSystem.Cartesian) {
                volumeDataset.dimX = dataSize[0];
                volumeDataset.dimY = dataSize[1];
                volumeDataset.dimZ = dataSize[2];
            } else {
                volumeDataset.dimX = gridSize[0];
                volumeDataset.dimY = gridSize[1];
                volumeDataset.dimZ = gridSize[2];
            }

            float[] densityData = processData();
            volumeDataset.data = densityData;

            Debug.Log("Loaded dataset in range: " + volumeDataset.GetMinDataValue() + "  -  " + volumeDataset.GetMaxDataValue());

            volumeDataset.FixDimensions();
            volumeDataset.rotation = Quaternion.Euler(90.0f, 0.0f, 0.0f);
        }

        private float[] processData() {

            float[] densities;

            // For grid-based Cartesian simulations, only need to load, flatten, and filter (if desired)
            if (simType == SimulationType.GridBased && coordinateSystem == CoordinateSystem.Cartesian) {

                float[,,] rawDensities = loadAndConvert3DData(dataset, dataSize[0], dataSize[1], dataSize[2]);
                densities = flatten3DData(rawDensities, dataSize[0], dataSize[1], dataSize[2]);
                if (filterData) {
                    densities = filter1DData(densities, filterValue);
                }
            // For the other types we need to bin the data, so those have all been grouped together
            } else {

                float[] r;
                float[] theta;
                float[] phi;
                float[] x;
                float[] y;
                float[] z;
                float[] densities1D;

                if (simType == SimulationType.GridBased) {

                    if (sphericalType == SphericalType.Uniform) {
                        // Generate r, theta, and phi arrays from given ranges
                        r = generateLerpArray(rMin, rMax, dataSize[0]);
                        theta = generateLerpArray(thetaMin, thetaMax, dataSize[1]);
                        phi = generateLerpArray(phiMin, phiMax, dataSize[2]);
                    } else {
                        // Just load the data in because we can't compute it
                        r = loadAndConvert1DData(rData, dataSize[0]);
                        theta = loadAndConvert1DData(thetaData, dataSize[1]);
                        phi = loadAndConvert1DData(phiData, dataSize[2]);
                    }

                    if (angleUnits == AngleUnits.Degrees) {
                        theta = convertToRadians(theta);
                        phi = convertToRadians(phi);
                    }

                    float[,,] rawDensities = loadAndConvert3DData(dataset, dataSize[0], dataSize[1], dataSize[2]);

                    // We flatten our key data arrays to unify with the 1D sources and because VolumeDataset
                    // needs the density as a 1D array anyway
                    x = new float[dataSize[0]*dataSize[1]*dataSize[2]];
                    y = new float[dataSize[0]*dataSize[1]*dataSize[2]];
                    z = new float[dataSize[0]*dataSize[1]*dataSize[2]];
                    densities1D = new float[dataSize[0]*dataSize[1]*dataSize[2]];

                    for (int i = 0; i < dataSize[0]; i++) {
                        for (int j = 0; j < dataSize[1]; j++) {
                            for (int k = 0; k < dataSize[2]; k++) {
                                x[i+j*dataSize[0]+k*dataSize[0]*dataSize[1]] = r[i]*Mathf.Sin(theta[j])*Mathf.Cos(phi[k]);
                                y[i+j*dataSize[0]+k*dataSize[0]*dataSize[1]] = r[i]*Mathf.Sin(theta[j])*Mathf.Sin(phi[k]);
                                z[i+j*dataSize[0]+k*dataSize[0]*dataSize[1]] = r[i]*Mathf.Cos(theta[j]);
                                densities1D[i+j*dataSize[0]+k*dataSize[0]*dataSize[1]] = rawDensities[i,j,k];
                            }
                        }
                    }
                } else {
                    if (coordinateSystem == CoordinateSystem.Cartesian) {

                        x = loadAndConvert1DData(xData, gridSize[0]*gridSize[1]*gridSize[2]);
                        y = loadAndConvert1DData(yData, gridSize[0]*gridSize[1]*gridSize[2]);
                        z = loadAndConvert1DData(zData, gridSize[0]*gridSize[1]*gridSize[2]);

                    } else {

                        r = loadAndConvert1DData(rData, gridSize[0]*gridSize[1]*gridSize[2]);
                        theta = loadAndConvert1DData(thetaData, gridSize[0]*gridSize[1]*gridSize[2]);
                        phi = loadAndConvert1DData(phiData, gridSize[0]*gridSize[1]*gridSize[2]);

                        x = new float[gridSize[0]*gridSize[1]*gridSize[2]];
                        y = new float[gridSize[0]*gridSize[1]*gridSize[2]];
                        z = new float[gridSize[0]*gridSize[1]*gridSize[2]];

                        for (int i = 0; i < gridSize[0]*gridSize[1]*gridSize[2]; i++) {
                            x[i] = r[i]*Mathf.Sin(theta[i])*Mathf.Cos(phi[i]);
                            y[i] = r[i]*Mathf.Sin(theta[i])*Mathf.Sin(phi[i]);
                            z[i] = r[i]*Mathf.Cos(theta[i]);
                        }
                    }
                    densities1D = loadAndConvert1DData(dataset, gridSize[0]*gridSize[1]*gridSize[2]);
                }
                densities = filterAndBinDensities(x, y, z, densities1D); 
            }
            return densities;
        }

        // Filters the density data (if desired) before binning it to fit in a Cartesian grid
        // Used on all supported data schemes except grid-based Cartesian!
        private float[] filterAndBinDensities(float[] x, float[] y, float[] z, float[] densitiesFlattened) {

            if (filterData) {
                ulong count = 0;
                for (int i = 0; i < densitiesFlattened.GetLength(0); i++) {
                    if (densitiesFlattened[i] >= filterValue) {
                        count++;
                    }
                }
                float[] xFiltered = new float[count];
                float[] yFiltered = new float[count];
                float[] zFiltered = new float[count];
                float[] densitiesFiltered = new float[count];

                Debug.Log(count);
                Debug.Log(densitiesFlattened.GetLength(0));

                ulong index = 0;
                for (int i = 0; i < densitiesFlattened.GetLength(0); i++) {
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

            float[] densitiesBinned = new float[gridSize[0]*gridSize[1]*gridSize[2]];
            int[] counts = new int[gridSize[0]*gridSize[1]*gridSize[2]];

            for (int i = 0 ; i < densitiesFlattened.GetLength(0); i++) {
                ushort xInd = (ushort)Mathf.Floor((x[i]-xMin)/cellX);
                ushort yInd = (ushort)Mathf.Floor((y[i]-yMin)/cellY);
                ushort zInd = (ushort)Mathf.Floor((z[i]-zMin)/cellZ);

                int dataInd = xInd+yInd*gridSize[0]+zInd*gridSize[0]*gridSize[1];

                densitiesBinned[dataInd] = (densitiesBinned[dataInd]*counts[dataInd] +
                                                    densitiesFlattened[i])/(counts[dataInd]+1);
                counts[dataInd]++;
            }

            return densitiesBinned;
        }

        // Returns a 1D vector of linearly interpolated values
        // (a la linspace in numpy, MATLAB, etc)
        private float[] generateLerpArray(float minValue, float maxValue, int numSteps) {
            float[] lerpArr = new float[numSteps];
            for (int i = 0; i < numSteps; i++) {
                    lerpArr[i] = Mathf.Lerp(minValue, maxValue, (float)i/(numSteps-1));
                }
            return lerpArr;
        }

        // Handles the different supported data types, pulls 1D data from HDF5 and casts to floats
        private float[] loadAndConvert1DData(string dataName, int dim) {
            float [] data = new float[dim];
            switch (contentFormat) 
                {
                    case DataContentFormat.Uint8:
                    {
                        byte[] temp = pull1DDataFromFile<byte>(dataName);
                        data = byteToFloat1DArray(temp, dim);
                        break;
                    }
                    case DataContentFormat.Uint16:
                    {
                        ushort[] temp = pull1DDataFromFile<ushort>(dataName);
                        data = ushortToFloat1DArray(temp, dim);
                        break;
                    }
                    case DataContentFormat.Uint32:
                    {
                        uint[] temp = pull1DDataFromFile<uint>(dataName);
                        data = uintToFloat1DArray(temp, dim);  
                        break;
                    }
                    case DataContentFormat.Int8:
                    {
                        sbyte[] temp = pull1DDataFromFile<sbyte>(dataName);
                        data = sbyteToFloat1DArray(temp, dim);
                        break;
                    }
                    case DataContentFormat.Int16:
                    {
                        short[] temp = pull1DDataFromFile<short>(dataName);
                        data = shortToFloat1DArray(temp, dim);
                        break;
                    }
                    case DataContentFormat.Int32:
                    {
                        int[] temp = pull1DDataFromFile<int>(dataName);
                        data = intToFloat1DArray(temp, dim);
                        break;
                    }
                    case DataContentFormat.Float32:
                    {
                        data = pull1DDataFromFile<float>(dataName);
                        break;
                    }
                    case DataContentFormat.Float64:
                    {
                        double[] temp = pull1DDataFromFile<double>(dataName);
                        data = doubleToFloat1DArray(temp, dim); 
                        break;                       
                    }
                    default: 
                    {
                        throw new NotImplementedException("Unimplemented data content format");
                    }
                }
            return data;
        }

        // Handles the different supported data types, pulls 3D data from HDF5 and casts to floats
        private float[,,] loadAndConvert3DData(string dataName, int firstDim, int secondDim, int thirdDim) {
            float [,,] data = new float[firstDim, secondDim, thirdDim];
            switch (contentFormat) 
                {
                    case DataContentFormat.Uint8:
                    {
                        byte[,,] temp = pull3DDataFromFile<byte>(dataName);
                        data = byteToFloat3DArray(temp, firstDim, secondDim, thirdDim);
                        break;
                    }
                    case DataContentFormat.Uint16:
                    {
                        ushort[,,] temp = pull3DDataFromFile<ushort>(dataName);
                        data = ushortToFloat3DArray(temp, firstDim, secondDim, thirdDim);
                        break;
                    }
                    case DataContentFormat.Uint32:
                    {
                        uint[,,] temp = pull3DDataFromFile<uint>(dataName);
                        data = uintToFloat3DArray(temp, firstDim, secondDim, thirdDim);  
                        break;
                    }
                    case DataContentFormat.Int8:
                    {
                        sbyte[,,] temp = pull3DDataFromFile<sbyte>(dataName);
                        data = sbyteToFloat3DArray(temp, firstDim, secondDim, thirdDim);
                        break;
                    }
                    case DataContentFormat.Int16:
                    {
                        short[,,] temp = pull3DDataFromFile<short>(dataName);
                        data = shortToFloat3DArray(temp, firstDim, secondDim, thirdDim);
                        break;
                    }
                    case DataContentFormat.Int32:
                    {
                        int[,,] temp = pull3DDataFromFile<int>(dataName);
                        data = intToFloat3DArray(temp, firstDim, secondDim, thirdDim);
                        break;
                    }
                    case DataContentFormat.Float32:
                    {
                        data = pull3DDataFromFile<float>(dataName);
                        break;
                    }
                    case DataContentFormat.Float64:
                    {
                        double[,,] temp = pull3DDataFromFile<double>(dataName);
                        data = doubleToFloat3DArray(temp, firstDim, secondDim, thirdDim);   
                        break;                     
                    }
                    default: 
                    {
                        throw new NotImplementedException("Unimplemented data content format");
                    }
                }
            return data;
        }

        // Uses HDF.PInvoke to set the data size and data type
        private void getDataPropertiesFromFile(string dataName) {
            // open the HDF5 file
            long file = H5F.open(filePath, H5F.ACC_RDONLY, H5P.DEFAULT);

            // open the relevant dataset
            long dset = H5D.open(file, dataName);

            // get the dataspace, pull the data dimensions from there
            long dspace = H5D.get_space(dset);
            int ndims = H5S.get_simple_extent_ndims(dspace);
            ulong[] dims = new ulong[ndims];
            H5S.get_simple_extent_dims(dspace, dims, null);

            // update the data size state
            for (int i = 0; i < ndims; i++) {
                dataSize[i] = (int)dims[i];
            }

            // get the datatype from the file
            long dtype = H5D.get_type(dset);
            long typecode = H5T.get_native_type(dtype, 0);

            // update the data type state
            contentFormat = GetFormatByTypeCode(typecode);

            H5T.close(dtype);
            H5S.close(dspace);
            H5D.close(dset);
            H5F.close(file);
        }

        private static DataContentFormat GetFormatByTypeCode(long typecode)
        {

            if (Convert.ToBoolean(H5T.equal(typecode, H5T.NATIVE_SHORT))) {
                return DataContentFormat.Int16;
            } else if (Convert.ToBoolean(H5T.equal(typecode, H5T.NATIVE_INT))) {
                return DataContentFormat.Int32;
            } else if (Convert.ToBoolean(H5T.equal(typecode, H5T.NATIVE_CHAR))) {
                return DataContentFormat.Int32;
            } else if (Convert.ToBoolean(H5T.equal(typecode, H5T.NATIVE_USHORT))) {
                return DataContentFormat.Int32;
            } else if (Convert.ToBoolean(H5T.equal(typecode, H5T.NATIVE_UINT))) {
                return DataContentFormat.Int32;
            } else if (Convert.ToBoolean(H5T.equal(typecode, H5T.NATIVE_UCHAR))) {
                return DataContentFormat.Int32;
            } else if (Convert.ToBoolean(H5T.equal(typecode, H5T.NATIVE_FLOAT))) {
                return DataContentFormat.Int32;
            } else if (Convert.ToBoolean(H5T.equal(typecode, H5T.NATIVE_DOUBLE))) {
                return DataContentFormat.Int32;
            } else {
                return DataContentFormat.Uint8;
            }
        }

        // Uses HDF.PInvoke to pull 1D data of a desired dimension from disk and into memory
        private T[] pull1DDataFromFile<T>(string dataName) {

            // open the HDF5 file
            long file = H5F.open(filePath, H5F.ACC_RDONLY, H5P.DEFAULT);

            // open the relevant dataset
            long dset = H5D.open(file, dataName);

            // get the datatype and associated number of bytes
            long dtype = H5D.get_type(dset);

            System.IntPtr size = H5T.get_size(dtype);

            // get the dataspace, pull the data dimensions from there
            long dspace = H5D.get_space(dset);
            int ndims = H5S.get_simple_extent_ndims(dspace);
            ulong[] dims = new ulong[ndims];
            H5S.get_simple_extent_dims(dspace, dims, null);

            T[] temp = new T[dims[0]];

            GCHandle pinnedArray = GCHandle.Alloc(temp, GCHandleType.Pinned);
            H5D.read(dset, dtype, H5S.ALL, H5S.ALL, H5P.DEFAULT, pinnedArray.AddrOfPinnedObject());
            pinnedArray.Free();

            H5T.close(dtype);
            H5S.close(dspace);
            H5D.close(dset);
            H5F.close(file);

            return temp;
        }

        // Uses HDF.PInvoke to pull 3D data of desired dimensions from disk and into memory
        private T[,,] pull3DDataFromFile<T>(string dataName) {

            // open the HDF5 file
            long file = H5F.open(filePath, H5F.ACC_RDONLY, H5P.DEFAULT);

            // open the relevant dataset
            long dset = H5D.open(file, dataName);

            // get the datatype and associated number of bytes
            long dtype = H5D.get_type(dset);

            System.IntPtr size = H5T.get_size(dtype);

            // get the dataspace, pull the data dimensions from there
            long dspace = H5D.get_space(dset);
            int ndims = H5S.get_simple_extent_ndims(dspace);
            ulong[] dims = new ulong[ndims];
            H5S.get_simple_extent_dims(dspace, dims, null);

            T[,,] temp = new T[dims[0], dims[1], dims[2]];

            GCHandle pinnedArray = GCHandle.Alloc(temp, GCHandleType.Pinned);
            H5D.read(dset, dtype, H5S.ALL, H5S.ALL, H5P.DEFAULT, pinnedArray.AddrOfPinnedObject());
            pinnedArray.Free();

            H5T.close(dtype);
            H5S.close(dspace);
            H5D.close(dset);
            H5F.close(file);

            return temp;
        }

        // Removes all values in a 1D array below a given value
        // Used for processing the grid-based Cartesian results, not for the others
        // because we flatten as part of the density binning process
        private float[] filter1DData(float[] data, float value) {
            ulong count = 0;
            for (int i = 0; i < data.GetLength(0); i++) {
                if (data[i] >= value) {
                    count++;
                }
            }
            float[] dataFiltered = new float[count];

            ulong index = 0;
            for (int i = 0; i < data.GetLength(0); i++) {
                if (data[i] >= value) {
                    dataFiltered[index] = data[i];
                    index++;
                }
            }
            return dataFiltered;
        }

        // Flattens a given 3D array into a 1D array in a predictable manner
        // Used for processing the grid-based Cartesian results, not for the others
        // because we flatten as part of the density binning process
        private float[] flatten3DData(float[,,] data, int firstDim, int secondDim, int thirdDim) {
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

        private float[] convertToRadians(float[] angles) {
            float[] radians = new float[angles.GetLength(0)];
            for (int i = 0; i < angles.GetLength(0); i++) {
                radians[i] = angles[i]*(Mathf.PI/180);
            }
            return radians;
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

        private float[] byteToFloat1DArray(byte[] data, int dim) {
            float[] converted = new float[dim];
            for (int i = 0; i < dim; i++) {
                converted[i] = (float)data[i];
            }
            return converted;
        }

        private float[] ushortToFloat1DArray(ushort[] data, int dim) {
            float[] converted = new float[dim];
            for (int i = 0; i < dim; i++) {
                converted[i] = (float)data[i];
            }
            return converted;
        }

        private float[] uintToFloat1DArray(uint[] data, int dim) {
            float[] converted = new float[dim];
            for (int i = 0; i < dim; i++) {
                converted[i] = (float)data[i];
            }
            return converted;
        }

        private float[] sbyteToFloat1DArray(sbyte[] data, int dim) {
            float[] converted = new float[dim];
            for (int i = 0; i < dim; i++) {
                converted[i] = (float)data[i];
            }
            return converted;
        }

        private float[] shortToFloat1DArray(short[] data, int dim) {
            float[] converted = new float[dim];
            for (int i = 0; i < dim; i++) {
                converted[i] = (float)data[i];
            }
            return converted;
        }

        private float[] intToFloat1DArray(int[] data, int dim) {
            float[] converted = new float[dim];
            for (int i = 0; i < dim; i++) {
                converted[i] = (float)data[i];
            }
            return converted;
        }

        private float[] doubleToFloat1DArray(double[] data, int dim) {
            float[] converted = new float[dim];
            for (int i = 0; i < dim; i++) {
                converted[i] = (float)data[i];
            }
            return converted;
        }

        private float[,,] byteToFloat3DArray(byte[,,] data, int firstDim, int secondDim, int thirdDim) {
            float[,,] converted = new float[firstDim, secondDim, thirdDim];
            for (int i = 0; i < firstDim; i++) {
                for (int j = 0; j < secondDim; j++) {
                    for (int k = 0; k < thirdDim; k++) {
                        converted[i,j,k] = (float)data[i,j,k];  
                    }
                }
            }
            return converted;
        }

        private float[,,] ushortToFloat3DArray(ushort[,,] data, int firstDim, int secondDim, int thirdDim) {
            float[,,] converted = new float[firstDim, secondDim, thirdDim];
            for (int i = 0; i < firstDim; i++) {
                for (int j = 0; j < secondDim; j++) {
                    for (int k = 0; k < thirdDim; k++) {
                        converted[i,j,k] = (float)data[i,j,k];  
                    }
                }
            }
            return converted;
        }

        private float[,,] uintToFloat3DArray(uint[,,] data, int firstDim, int secondDim, int thirdDim) {
            float[,,] converted = new float[firstDim, secondDim, thirdDim];
            for (int i = 0; i < firstDim; i++) {
                for (int j = 0; j < secondDim; j++) {
                    for (int k = 0; k < thirdDim; k++) {
                        converted[i,j,k] = (float)data[i,j,k];  
                    }
                }
            }
            return converted;
        }

        private float[,,] sbyteToFloat3DArray(sbyte[,,] data, int firstDim, int secondDim, int thirdDim) {
            float[,,] converted = new float[firstDim, secondDim, thirdDim];
            for (int i = 0; i < firstDim; i++) {
                for (int j = 0; j < secondDim; j++) {
                    for (int k = 0; k < thirdDim; k++) {
                        converted[i,j,k] = (float)data[i,j,k];  
                    }
                }
            }
            return converted;
        }

        private float[,,] shortToFloat3DArray(short[,,] data, int firstDim, int secondDim, int thirdDim) {
            float[,,] converted = new float[firstDim, secondDim, thirdDim];
            for (int i = 0; i < firstDim; i++) {
                for (int j = 0; j < secondDim; j++) {
                    for (int k = 0; k < thirdDim; k++) {
                        converted[i,j,k] = (float)data[i,j,k];  
                    }
                }
            }
            return converted;
        }

        private float[,,] intToFloat3DArray(int[,,] data, int firstDim, int secondDim, int thirdDim) {
            float[,,] converted = new float[firstDim, secondDim, thirdDim];
            for (int i = 0; i < firstDim; i++) {
                for (int j = 0; j < secondDim; j++) {
                    for (int k = 0; k < thirdDim; k++) {
                        converted[i,j,k] = (float)data[i,j,k];  
                    }
                }
            }
            return converted;
        }

        private float[,,] doubleToFloat3DArray(double[,,] data, int firstDim, int secondDim, int thirdDim) {
            float[,,] converted = new float[firstDim, secondDim, thirdDim];
            for (int i = 0; i < firstDim; i++) {
                for (int j = 0; j < secondDim; j++) {
                    for (int k = 0; k < thirdDim; k++) {
                        converted[i,j,k] = (float)data[i,j,k];  
                    }
                }
            }
            return converted;
        }
    }
}
