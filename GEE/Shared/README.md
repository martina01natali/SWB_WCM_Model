# README file for GEE - Shared folder
# Notebooks for download of $\sigma^0$ and NDVI values from Sentinel-1, Sentinel-2 GEE collections

## Dependencies 

These codes require the installation of the Earth Engine API, `ee`. You can find more info on the installation procedure here: [Python installation of GEE](https://developers.google.com/earth-engine/guides/python_install). \
This code runs on browser-based notebooks only (Google Colaboratory, Jupyter Notebooks, etc...). \
Be aware that you won't need to install the Google Cloud APK to run the code.

The full dependencies required are provided in the file of the environment google, `google.yml`. To install the virtual environment by using conda, run the lines below in the terminal after having activated `conda`.

```python
conda env create -f environment.yml # install environment from file environment.yml
conda activate <myenv> # activate environment with name provided by header of .yml file
conda info --envs # to check if everything worked fine
```
If you don't use `conda`, manually install the dependencies listed in the environment file by using `pip`.