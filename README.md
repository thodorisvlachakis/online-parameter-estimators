# Simulation and Modeling – Assignment 2
**Topic:** Real-Time Parameter Estimation (Gradient Descent Method (Online Version), Lyapunov Method)
**Course:** Simulation and Modeling of Dynamic Systems  

## 📘 Description
This project focuses on **online estimation of unknown parameters** for dynamic systems using **real-time methods**.

The project focuses on:
1. **Gradient Descent (online) Method** – an adaptive estimator that updates parameters in real-time based on the observed system response.  
2. **Lyapunov Method** – a method ensuring stability of parameter estimates while adapting in real-time.  

Two systems are studied:

1. **Mass-Spring-Damper System (MSD)**  
   The system is described by:
   
   $$ m\ddot{x}(t) + b\dot{x}(t) + kx(t) = u(t) $$

   where:
   - \(x(t)\) [m]: displacement  
   - \(m > 0\): mass  
   - \(b > 0\): damping coefficient  
   - \(k > 0\): spring constant  
   - \(u(t)\): external force  

   The goal is to estimate the unknown parameters \(m\), \(b\), and \(k\) in real-time using:  
   - **Gradient (online) method** for different inputs \(u(t)\)  
   - **Lyapunov method** with parallel and mixed structures  
   - Analysis of estimation performance under **output noise**  

3. **Nonlinear Roll Angle System**  
   The system describes the roll dynamics of an aircraft:

   $$  \ddot{r}(t) = -a_1 \dot{r}(t) - a_2 \sin(r(t)) + a_3 \dot{r}^2(t) \sin(2r(t)) + b u(t) + d(t)  $$

   where:
   - \(r(t)\) [rad]: roll angle  
   - \(a_i > 0, b > 0\): unknown system parameters  
   - \(u(t)\): control input  
   - \(d(t)\): external disturbances  

   The objectives are:  
   - Implement a **feedback controller** to track a smooth reference trajectory \(r_d(t)\)  
   - Apply a **Lyapunov-based online estimator** for the unknown parameters  
   - Study the effect of external disturbances on estimation accuracy  

---

## ⚙️ Implementation
- **Language:** MATLAB  
- **Approach:**
  1. Define system models and parameter estimation logic in **functions** (`functions/` folder)
  2. Implement **both Gradient Descent and Lyapunov methods** for real-time parameter estimation
  3. Execute each exercise script (`exercises/` folder) to simulate the system, run estimators, and generate plots  
  4. Generate plots of:
     - System response vs. estimated response  
     - Parameter estimates over time  
     - Estimation error
  5. Compare performance of **Gradient Descent** vs **Lyapunov** estimators  
  6. Study effects of noise and disturbances where applicable  

---

## 📂 Repository Structure
code/
├── functions/
│ ├── composite_stricture_model_adaptive_estimation_MSD.m
│ ├── composite_stricture_model_adaptive_estimation_RollAngleSys.m
│ ├── desiredRollAngleTrajectory.m
│ ├── gradient_descent_online_estimator_MSD.m
│ ├── parallel_structure_model_adaptive_estimation_MSD.m
│ ├── rollAngleDynamicSystem.m
│ ├── rollAngleSysController.m
│ └── systemEquationsOfState.m
└── exercises/
├── Exercise_1a_Parameter_Estimation_MassSpringDamper_System.m
├── Exercise_1bc_Parameter_Estimation_MassSpringDamper_System.m
├── Exercise_2a_Roll_Angle_Sys_Controller_Implementation.m
└── Exercise_2bc_Parameter_Estimation_Roll_Angle_Sys.m

statement/
report/
figures_and_results/
tests/


- `functions/` contains all MATLAB functions implementing system models, controllers, and online estimators  
- `exercises/` contains scripts for each exercise that run simulations and produce plots/results
- `statement/` contains the assignment description PDF file
- `report/` contains the report PDF file with mathematical analyses, comments and conclusions
- `figures_and_results/` stores generated plot and screenshots for reference or report inclusion  
- `tests/` is optional for temporary scripts and experiments  

---

## ▶️ How to Run
1. Open MATLAB and navigate to the `code/` folder.  
2. Run any exercise script from `exercises/` to simulate the corresponding system and estimator.  
3. The scripts will generate plots showing:  
   - System response vs. estimated response  
   - Estimation errors  
   - Parameter evolution over time  
4. Optionally, save plots to `figures_and_results/` for reporting purposes.  

---

## 📊 Results
- The figures demonstrate the performance of **online parameter estimators** for different inputs and system configurations  
- Comparisons between **gradient method** and **Lyapunov method**, as well as the effect of noise/disturbances, are included in the plots  
- Estimated parameters should converge to their true values, and estimation errors can be analyzed for accuracy assessment  
