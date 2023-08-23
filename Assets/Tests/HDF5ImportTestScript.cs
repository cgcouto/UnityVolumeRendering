using System.Collections;
using System.Collections.Generic;
using NUnit.Framework;
using UnityEngine;
using UnityEngine.TestTools;
// using UnityVolumeRendering;

namespace Tests
{
    public class HDF5ImportTestScript
    {
        // A Test behaves as an ordinary method
        [Test]
        public void NewTestScriptSimplePasses()
        {
            // Use the Assert class to test conditions
            // float rMin = 0.0f;
            // float rMax = 10.0f;
            // float thetaMin = 0;
            // float thetaMax = 180;
            // float phiMin = 0;
            // float phiMax = 360;
            // DataContentFormat format = DataContentFormat.Float64;
            // CoordinateSystem system = CoordinateSystem.Spherical;
            // AngleUnits angles = AngleUnits.Degrees;
            // int dimX = 10;
            // int dimY = 10;
            // int dimZ = 10;
            // int gridX = 20;
            // int gridY = 20;
            // int gridZ = 20;
            // double[,,] data = new double[dimX,dimY,dimZ];
            // for (int j = 0; j < 10; j++) {
            //     for (int k = 0; k < 10; k++) {
            //         data[5,j,k] = 100.0;
            //     }
            // }
            // HDF5DatasetImporter importer = new HDF5DatasetImporter("", "", dimX, dimY, dimZ, format,  
            //                                                         system, angles, rMin, rMax, thetaMin, thetaMax,
            //                                                         phiMin, phiMax, gridX, gridY, gridZ);
            // VolumeDataset volumeDataset = importer.ImportAsync();


        }

        // A UnityTest behaves like a coroutine in Play Mode. In Edit Mode you can use
        // `yield return null;` to skip a frame.
        [UnityTest]
        public IEnumerator NewTestScriptWithEnumeratorPasses()
        {
            // Use the Assert class to test conditions.
            // Use yield to skip a frame.
            yield return null;
        }
    }
}
