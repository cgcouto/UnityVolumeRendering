using UnityEngine;
using UnityEditor;
using System.IO;
using System;
using System.Threading.Tasks;

namespace UnityVolumeRendering
{
    /// <summary>
    /// Editor window for importing datasets.
    /// </summary>
    public class HDF5DatasetImporterEditorWindow : EditorWindow
    {
        private string fileToImport;

        private string dataset;

        private int dimX;
        private int dimY;
        private int dimZ;

        private float rMin;
        private float rMax;
        private float thetaMin;
        private float thetaMax;
        private float phiMin;
        private float phiMax;

        private int gridX;
        private int gridY;
        private int gridZ;

        private bool filterToggle;
        private float filterLessThan;

        private DataContentFormat dataFormat = DataContentFormat.Uint32;

        private bool importing = false;

        private CoordinateSystem coordinateSystem = CoordinateSystem.Cartesian;

        private AngleUnits angleUnits = AngleUnits.Degrees;


        public void Initialise(string filePath)
        {
            fileToImport = filePath;

            if (Path.GetExtension(fileToImport) == ".ini")
                fileToImport = fileToImport.Substring(0, fileToImport.Length - 4);

            // Try parse ini file (if available)
            DatasetIniData iniData = DatasetIniReader.ParseIniFile(fileToImport + ".ini");
            if (iniData != null)
            {
                dimX = iniData.dimX;
                dimY = iniData.dimY;
                dimZ = iniData.dimZ;
                dataFormat = iniData.format;
                dataset = iniData.dataset;
                rMin = iniData.rMin;
                rMax = iniData.rMax;
                thetaMin = iniData.thetaMin;
                thetaMax = iniData.thetaMax;
                phiMin = iniData.phiMin;
                phiMax = iniData.phiMax;
                gridX = iniData.gridX;
                gridY = iniData.gridY;
                gridZ = iniData.gridZ;
                filterLessThan = iniData.filterLessThan;
            }
            // this.position.width = 400.0f;
            this.minSize = new Vector2(300.0f, 200.0f);
        }

        private async Task ImportDatasetAsync()
        {
            using (ProgressHandler progressHandler = new ProgressHandler(new EditorProgressView(), "HDF5 import"))
            {
                progressHandler.ReportProgress(0.0f, "Importing HDF5 dataset");

                HDF5DatasetImporter importer = new HDF5DatasetImporter();

                int[] dataSize = {dimX, dimY, dimZ};

                if (coordinateSystem == CoordinateSystem.Cartesian) {
                    importer = new HDF5DatasetImporter(fileToImport, dataset, dataSize, dataFormat, coordinateSystem);
                } else {
                    int[] gridSize = {gridX, gridY, gridZ};
                    importer = new HDF5DatasetImporter(fileToImport, dataset, dataSize, dataFormat,  
                                                                            coordinateSystem, angleUnits, rMin, rMax, thetaMin, thetaMax,
                                                                            phiMin, phiMax, gridSize, filterToggle, filterLessThan);
                }

                VolumeDataset volumeDataset = await importer.ImportAsync();

                if (dataset != null)
                {
                    if (EditorPrefs.GetBool("DownscaleDatasetPrompt"))
                    {
                        if (EditorUtility.DisplayDialog("Optional DownScaling",
                            $"Do you want to downscale the dataset? The dataset's dimension is: {volumeDataset.dimX} x {volumeDataset.dimY} x {volumeDataset.dimZ}", "Yes", "No"))
                        {
                            Debug.Log("Async dataset downscale. Hold on.");
                            progressHandler.ReportProgress(0.7f, "Downscaling dataset");
                            await Task.Run(() =>  volumeDataset.DownScaleData());
                        }
                    }
                    progressHandler.ReportProgress(0.8f, "Creating object");
                    VolumeRenderedObject obj = await VolumeObjectFactory.CreateObjectAsync(volumeDataset);
                }
                else
                {
                    Debug.LogError("Failed to import datset");
                }

                this.Close();
            }
        }

        private async void StartImport()
        {
            try
            {
                importing = true;
                await ImportDatasetAsync();
            }
            catch (Exception ex)
            {
                importing = false;
                Debug.LogException(ex);
            }
            importing = false;
        }

        private void OnGUI()
        {
            if (importing)
            {
                EditorGUILayout.LabelField("Importing dataset. Please wait..");
            }
            else
            {
                EditorGUIUtility.labelWidth = 0.5f*this.position.width;
                dataset = EditorGUILayout.TextField("Dataset", dataset);
                dimX = EditorGUILayout.IntField("First dimension of current data (X or r)", dimX);
                dimY = EditorGUILayout.IntField("Second dimension of current data (Y or θ)", dimY);
                dimZ = EditorGUILayout.IntField("Third dimension of current data (Z or Φ)", dimZ);
                dataFormat = (DataContentFormat)EditorGUILayout.EnumPopup("Data format", dataFormat);
                coordinateSystem = (CoordinateSystem)EditorGUILayout.EnumPopup("Coordinate system", coordinateSystem);

                if (coordinateSystem == CoordinateSystem.Spherical) {
                    rMin = EditorGUILayout.FloatField("Minimum r value", rMin);
                    rMax = EditorGUILayout.FloatField("Maximum r value", rMax);
                    thetaMin = EditorGUILayout.FloatField("Minimum θ value", thetaMin);
                    thetaMax = EditorGUILayout.FloatField("Maximum θ value", thetaMax);
                    phiMin = EditorGUILayout.FloatField("Minimum Φ value", phiMin);
                    phiMax = EditorGUILayout.FloatField("Maximum Φ value", phiMax);
                    angleUnits = (AngleUnits)EditorGUILayout.EnumPopup("Angle units", angleUnits);
                    gridX = EditorGUILayout.IntField("X dimension of Cartesian grid", gridX);
                    gridY = EditorGUILayout.IntField("Y dimension of Cartesian grid", gridY);
                    gridZ = EditorGUILayout.IntField("Z dimension of Cartesian grid", gridZ); 
                    filterToggle = EditorGUILayout.Toggle("Filter out densities less than some value?", filterToggle);

                    if (filterToggle) {
                        filterLessThan = EditorGUILayout.FloatField("Filter value", filterLessThan);
                    }
                    
                }

                if (GUILayout.Button("Import"))
                {
                    StartImport();
                }

                if (GUILayout.Button("Cancel"))
                    this.Close();
            }
        }
    }
}
