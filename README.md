# Simulating WGS cfDNA with scWGS
A python package that simulates cfDNA data by mixing scWGS data. 

## Getting started 

1. Create a working directory:
   ```
   mkdir -p path/to/project
   cd path/to/project
   ```
2. Clone the repository through git:
      ```
      git clone --depth 1 https://github.com/matteolepur/cfsim.git
      ```
 
3. Download [pixi](https://pixi.prefix.dev/latest/), version >0.72.0:
    ```
    pixi init
    ```
4. You will be prompted to extended the existing `pyproject.toml` click yes:
   ```
   A 'pyproject.toml' file already exists. Do you want to extend it with the '[tool.pixi]' configuration? [y/N]
   y
   ```
4. Verify the installation worked:
    ```
    pixi run cfsim
    ```
--------
### Commands
- `cfsim simulate`

### Inputs 

### Outputs

#### Last time
- implement plotting code
