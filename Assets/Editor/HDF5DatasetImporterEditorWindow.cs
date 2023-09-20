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

        private string rData;
        private string thetaData;
        private string phiData;

        private string xData;
        private string yData;
        private string zData;

        private SphericalType sphericalType;

        private bool importing = false;

        private CoordinateSystem coordinateSystem = CoordinateSystem.Cartesian;

        private AngleUnits angleUnits = AngleUnits.Degrees;

        private SimulationType simType = SimulationType.GridBased;


        public void Initialise(string filePath)
        {
            fileToImport = filePath;

            if (Path.GetExtension(fileToImport) == ".ini")
                fileToImport = fileToImport.Substring(0, fileToImport.Length - 4);

            // Try parse ini file (if available)
            DatasetIniData iniData = DatasetIniReader.ParseIniFile(fileToImport + ".ini");
            if (iniData != null)
            {
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
                filterToggle = iniData.filterBool;
                rData = iniData.rData;
                thetaData = iniData.thetaData;
                phiData = iniData.phiData;
                xData = iniData.xData;
                yData = iniData.yData;
                zData = iniData.zData;
                coordinateSystem = iniData.coordinateSystem;
                simType = iniData.simType;
                sphericalType = iniData.sphericalType;
                angleUnits = iniData.angleUnits;
            }
            this.minSize = new Vector2(300.0f, 200.0f);
        }

        private async Task ImportDatasetAsync()
        {
            using (ProgressHandler progressHandler = new ProgressHandler(new EditorProgressView(), "HDF5 import"))
            {
                progressHandler.ReportProgress(0.0f, "Importing HDF5 dataset");

                HDF5DatasetImporter importer = new HDF5DatasetImporter();

                int[] gridSize = {gridX, gridY, gridZ};

                if (coordinateSystem == CoordinateSystem.Cartesian && simType == SimulationType.ParticleBased) {
                    importer = new HDF5DatasetImporter(fileToImport, dataset, simType, coordinateSystem,
                                                        angleUnits, xData, yData, zData, gridSize, filterToggle, filterLessThan);
                } else if (coordinateSystem == CoordinateSystem.Spherical && simType == SimulationType.GridBased 
                           && sphericalType == SphericalType.Uniform) {
                    importer = new HDF5DatasetImporter(fileToImport, dataset, simType, coordinateSystem, 
                                                        sphericalType, angleUnits, rMin, rMax, thetaMin, thetaMax,
                                                        phiMin, phiMax, gridSize, filterToggle, filterLessThan);
                } else if (coordinateSystem == CoordinateSystem.Spherical && simType == SimulationType.GridBased
                           && sphericalType == SphericalType.NonUniform) {
                    importer = new HDF5DatasetImporter(fileToImport, dataset, simType, coordinateSystem, 
                                                        sphericalType, angleUnits, rData, thetaData, phiData,
                                                        gridSize, filterToggle, filterLessThan);                    
                } else if (coordinateSystem == CoordinateSystem.Spherical && simType == SimulationType.ParticleBased) {
                    importer = new HDF5DatasetImporter(fileToImport, dataset, simType, coordinateSystem,
                                                        rData, thetaData, phiData, angleUnits, gridSize,filterToggle, filterLessThan);
                } else {
                    importer = new HDF5DatasetImporter(fileToImport, dataset, simType, coordinateSystem,
                                                       filterToggle,filterLessThan);
                }

                VolumeDataset volumeDataset = await importer.ImportAsync();

                if (dataset != null)
                {
                    if (EditorPrefs.GetBool("DownscaleDatasetPrompt"))
                    {
                        if (EditorUtility.DisplayDialog("Optional DownScaling",
                            $"Do you want to downscale the dataset?", "Yes", "No"))
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
                DrawUILine(Color.gray);
                simType = (SimulationType)EditorGUILayout.EnumPopup("Simulation type", simType);
                if (simType == SimulationType.GridBased) {
                    dataset = EditorGUILayout.TextField("Dataset for density data", dataset);
                } else if (simType == SimulationType.ParticleBased) {
                    dataset = EditorGUILayout.TextField("Dataset for density data", dataset);
                    if (coordinateSystem == CoordinateSystem.Cartesian) {
                        xData = EditorGUILayout.TextField("Dataset for positions in X", xData);
                        yData = EditorGUILayout.TextField("Dataset for positions in Y", yData);
                        zData = EditorGUILayout.TextField("Dataset for positions in Z", zData);
                    } else if (coordinateSystem == CoordinateSystem.Spherical) {
                        rData = EditorGUILayout.TextField("Dataset for positions in r", rData);
                        thetaData = EditorGUILayout.TextField("Dataset for positions in θ", thetaData);
                        phiData = EditorGUILayout.TextField("Dataset for positions in Φ", phiData);                        
                    }
                }
                DrawUILine(Color.gray);
                coordinateSystem = (CoordinateSystem)EditorGUILayout.EnumPopup("Coordinate system", coordinateSystem);

                if (coordinateSystem == CoordinateSystem.Spherical) {
                    if (simType == SimulationType.GridBased) {
                        sphericalType = (SphericalType)EditorGUILayout.EnumPopup("Type of spherical grid", sphericalType);
                        if (sphericalType == SphericalType.Uniform) {
                            EditorGUILayout.SelectableLabel("This option works when your r, θ, and Φ values are evenly spaced!");
                            rMin = EditorGUILayout.FloatField("Minimum r value", rMin);
                            rMax = EditorGUILayout.FloatField("Maximum r value", rMax);
                            thetaMin = EditorGUILayout.FloatField("Minimum θ value", thetaMin);
                            thetaMax = EditorGUILayout.FloatField("Maximum θ value", thetaMax);
                            phiMin = EditorGUILayout.FloatField("Minimum Φ value", phiMin);
                            phiMax = EditorGUILayout.FloatField("Maximum Φ value", phiMax);
                        } else if (sphericalType == SphericalType.NonUniform) {
                            EditorGUILayout.SelectableLabel("If your data has non-uniform spacing going on, use this!");
                            rData = EditorGUILayout.TextField("Dataset for r values", rData);
                            thetaData = EditorGUILayout.TextField("Dataset for θ values", thetaData);
                            phiData = EditorGUILayout.TextField("Dataset for Φ values", phiData);
                        }
                    }
                    angleUnits = (AngleUnits)EditorGUILayout.EnumPopup("Angle units", angleUnits);
                }
                DrawUILine(Color.gray);

                if (coordinateSystem == CoordinateSystem.Spherical || simType == SimulationType.ParticleBased) {
                    gridX = EditorGUILayout.IntField("X dimension of desired Cartesian grid", gridX);
                    gridY = EditorGUILayout.IntField("Y dimension of desired Cartesian grid", gridY);
                    gridZ = EditorGUILayout.IntField("Z dimension of desired Cartesian grid", gridZ); 
                }

                filterToggle = EditorGUILayout.Toggle("Filter out densities less than some value?", filterToggle);

                if (filterToggle) {
                    filterLessThan = EditorGUILayout.FloatField("Filter value", filterLessThan);
                }  
                DrawUILine(Color.gray);

                if (GUILayout.Button("Import"))
                {
                    StartImport();
                }

                if (GUILayout.Button("Cancel"))
                    this.Close();
            }
        }
        // from post #6 on https://forum.unity.com/threads/horizontal-line-in-editor-window.520812/
        private static void DrawUILine(Color color, int thickness = 2, int padding = 10)
        {
            Rect r = EditorGUILayout.GetControlRect(GUILayout.Height(padding+thickness));
            r.height = thickness;
            r.y+=padding/2;
            r.x-=2;
            r.width +=6;
            EditorGUI.DrawRect(r, color);
        }
    }
}
