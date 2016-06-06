# isingmodel

To output all matrix:
  - Compile file:
  ```
  $ make picture=1
  ```
  - Output data:
    - For 2d ising model:
    ```
    $ ./2dising.sh
    ```
    (If it doesn't work, type in:
    ```
    $ chmod +x ./2dising.sh
    ```
    ).
    - For 3d ising model:
    ```
    $ ./3dising.sh
    ```
    (If it doesn't work, type in:
    ```
    $ chmod +x ./3dising.sh
    ```
    ).

To output energy, heat capacity, magnetization, magnetization susceptibility:
  - Compile file:
  ```
  $ make
  ```
  - Output data:
    - For 2d ising model:
    ```
    $ ./2dising.sh
    ```
    - For 3d ising model:
    ```
    $ ./3dising.sh
    ```

- 2d ising model data will be in ./data/2d
- 3d ising model data will be in ./data/3d

To clean all data:
```sh
$ make cleandata
```
To clean 2d data:
```sh
$ make clean2ddata
```
To clean 3d data:
```sh
$ make clean3ddata
```

How to access Titan box:
```sh
$ ssh username@ieng6.ucsd.edu
$ ssh username@igpu6-210.ucsd.edu
```

How to run the cuda code:
  - Compile file:
  ```
  $ make picture=1
  ```
  - Output data:
  ```
  $ ./2dising.sh
  ```
  (If it doesn't work, type in:
  ```
  $ chmod +x ./2dising_cuda.sh
  ```
  ).
