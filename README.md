# RRbot-Robust-Control

## Aim
The aim is to design a Robust Control law for the RRbot and perform trajectory tracking.

## Design Procedure
The desired cubic polynomial trajectories are calculated based on the state parameters given for the initial and final conditions.

The control law is designed in MATLAB, and the parameters are tuned to achieve desirable state trajectories and control input trajectories. Once the results are satisfactory, further experiments are performed in Gazebo.

## Results

- The implemented control law performed as expected to perform trajectory tracking within the given limited torque constraints. 

- The state and control input trajectories without a boundary layer are captured as follows. Though, the robot converges to the desired trajectories, the chattering issue in the control inputs can be observed. 

<img width="740" alt="Screenshot 2023-05-27 at 11 27 06 PM" src="https://github.com/kt-krutarthtrivedi/RRbot-Robust-Control/assets/134632027/706e6144-8e3c-4748-883a-68f7f952a0be">

&nbsp; 

- The state and control input trajectories with a boundary layer are captured as follows. As seen, the chattering in the control inputs is completely removed and we have desirable control input trajectories.

<img width="740" alt="Screenshot 2023-05-27 at 11 34 21 PM" src="https://github.com/kt-krutarthtrivedi/RRbot-Robust-Control/assets/134632027/267ca53c-abde-4f2c-9214-f67ce6eaa3e9">


## Report
Please find a detailed description of the control law design and results in the [report](https://github.com/kt-krutarthtrivedi/RRbot-Robust-Control/blob/main/media/Report.pdf).
