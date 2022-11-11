# CourseAlyaACC

Repo for the course "AlyaACC" at BSC-CASE. Contains the following:

- `lectureAlyaGPU.pdf`: slides for the course;
- `src/`: source code for the course;
  - `vecAdd`: vector addition example;
  - `conjGrad/`: source code for the conjugate gradient method;
  - `FEM1D/`: source code for the 1D FEM example;

Instructions on how to access CTE-POWER and copy this repo to your home directory can be found on slide 2 of `lectureAlyaGPU.pdf`.

## Building the examples

The examples are built using Make. To build the examples, run the following commands:

```bash
cd src/[example]
make
```

Makefiles are provided for each example, and can be viewed for insight into how the examples are built.