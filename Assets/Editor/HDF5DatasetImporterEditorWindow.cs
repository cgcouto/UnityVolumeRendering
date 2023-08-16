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

        private DataContentFormat dataFormat = DataContentFormat.Uint32;

        private bool importing = false;

        public void Initialise(string filePath)
        {
            fileToImport = filePath;

            if (Path.GetExtension(fileToImport) == ".ini")
                fileToImport = fileToImport.Substring(0, fileToImport.Length - 4);

            // // Try parse ini file (if available)
            // DatasetIniData initData = DatasetIniReader.ParseIniFile(fileToImport + ".ini");
            // if (initData != null)
            // {
            //     dimX = initData.dimX;
            //     dimY = initData.dimY;
            //     dimZ = initData.dimZ;
            //     bytesToSkip = initData.bytesToSkip;
            //     dataFormat = initData.format;
            //     endianness = initData.endianness;
            // }

            this.minSize = new Vector2(300.0f, 200.0f);
        }

        private async Task ImportDatasetAsync()
        {
            using (ProgressHandler progressHandler = new ProgressHandler(new EditorProgressView(), "HDF5 import"))
            {
                progressHandler.ReportProgress(0.0f, "Importing HDF5 dataset");

                HDF5DatasetImporter importer = new HDF5DatasetImporter(fileToImport, dataset, dimX, dimY, dimZ, dataFormat);
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
                dataset = EditorGUILayout.TextField("Dataset", dataset);
                dimX = EditorGUILayout.IntField("X dimension", dimX);
                dimY = EditorGUILayout.IntField("Y dimension", dimY);
                dimZ = EditorGUILayout.IntField("Z dimension", dimZ);
                dataFormat = (DataContentFormat)EditorGUILayout.EnumPopup("Data format", dataFormat);

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
