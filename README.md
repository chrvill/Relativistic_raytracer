(WIP)

# Relativistic-renderer
![Kerr black hole](./renders/kerr.png)

C++ code for performing ray tracing through curved spacetime. The focus is on rendering rotating black holes, described by the Kerr metric.

## The physics

Ray tracing through curved spacetime means tracing geodesics. These are the paths through spacetime followed by any free-falling particle, and they maximize the time that an observer following the path measures.

Geodesics are described usign the [geodesic equation](https://en.wikipedia.org/wiki/Geodesics_in_general_relativity):

$$
\begin{equation}
    \frac{d^2 x^\mu}{d\lambda^2} + \Gamma^\mu_{\rho \sigma} \frac{dx^\rho}{d\lambda} \frac{dx^\sigma}{d\lambda} = 0
\end{equation}
$$

for $\mu = 0, 1, 2, 3$. We use the [Einstein notation convention](https://en.wikipedia.org/wiki/Einstein_notation), meaning that we implicitly sum over $\rho$ and $\sigma$. Here $\lambda$ is an affine parameter for the geodesic, while $\Gamma^\mu_{\rho \sigma}$ are the [Christoffel symbols](https://en.wikipedia.org/wiki/Christoffel_symbols), given by 

$$
\begin{equation}
    \Gamma^\mu_{\rho \sigma} = \frac{1}{2}g^{\mu \nu}\left(\partial_\rho g_{\sigma \nu} + \partial_\sigma g_{\rho \nu} - \partial_\nu g_{\rho \sigma}\right)
\end{equation}
$$

and $g_{\mu \nu}$ is the [metric tensor](https://en.wikipedia.org/wiki/Metric_tensor), while $g^{\mu \nu}$ is its inverse. In (3 + 1)-dimensional spacetime all the indices take on four values, so the second term on the left hand side of (1) in reality represents $4^2 = 16$ different terms. We have to calculate all of these Christoffel symbols for the metric we are considering in order to trace out geodesics in the spacetime. But luckily, the Christoffel symbols are symmetric in the lower indices, meaning that $\Gamma^\mu_{\rho \sigma} = \Gamma^\mu_{\sigma \rho}$. This reduces the total number of distinct Christoffel symbols by half, but the procedure can in general be very tedious to perform. Here we use the code provided [here](https://github.com/chrvill/Black-hole-GR) to derive expressions for all the Christoffel symbols for a given metric.

#### Camera rays

We generate rays from the camera in a standard way, but the important thing to note is that these rays are emitted in the rest frame of the camera. Meanwhile, the geodesic tracing requires us to work in the global frame describing the spacetime. So we have to convert from the rest frame of the camera to the global frame. 

Assume the camera has four-velocity $u^\mu$ measured in the global frame. Then its four-velocity measured in some local (inertial) frame at the position of the camera is given by 

$$
\begin{equation}
    \tilde{u}^m = e^{m}_{\;\;\mu} u^\mu
\end{equation}
$$

where $e^{m}_{\;\;\mu}$ are [tetrads/frame fields](https://en.wikipedia.org/wiki/Frame_fields_in_general_relativity) describing the local frame (which is to say, each observer in relative motion to each other will have a different set of frame fields). Since $\tilde{u}^m$ is a four-velocity in an inertial frame it can also be written as $\tilde{u}^m = \gamma\left(1, \mathbf{v}\right)$, where $\gamma = \frac{1}{\sqrt{1 - \mathbf{v}^2}}$, and $\mathbf{v}$ is the velocity of the camera in the local inertial frame. In order to transform the camera rays from the rest frame of the camera to the local inertial frame specified by $e^{m}_{\;\;\mu}$ we can then simply Lorentz transform them:

$$
\begin{equation}
    p'^m = \Lambda^{m}_{\;\;n}(\mathbf{v})p^n
\end{equation}
$$

where $\Lambda^m_{\;\;n}(\mathbf{v})$ is this Lorentz transformation, $p^n$ is the four-momentum of a camera ray in the camera rest frame, and $p'^m$ is the four-momentum of the camera ray in the other local inertial frame. 

Then finally we use the inverse tetrads $e^{*\;\mu}_m$ to transform from the local inertial frame to the global frame:

$$
\begin{equation}
    p^\mu = e^{*\;\mu}_m p'^m
\end{equation}
$$

where $p^\mu$ is the four-momentum of a camera ray as measured in the global frame. The reason for going through the intermediate step of transforming to another local inertial frame, rather than just transforming directly from the camera's rest frame to the global frame, is that we typically choose some clever inertial frame where the tetrads take particularly nice forms. In Schwarzschild, for example, the choice is typically the standard shell observers, while in Kerr we typically choose Zero Angular Momentum Observers (ZAMOs).  

#### The numerical integration scheme

In $n$-dimensional spacetime the geodesic equation produces $n$ coupled, second order, ordinary differential equations. Solving these equations is not any different from solving any other set of coupled differential equations - one can apply a standard numerical integration scheme. An important consideration here, though, is that we are going to be running into coordinate singularities when nearing the event horizon of the black holes we render. So the numerical scheme one chooses has to be robust enough to properly handle those divergences. Working in spherical/Boyer-Lindquist coordinates, we will also encounter a coordinate singularity $\theta = 0$ at the vertical axis, which will produce similar effects. 

Here I've used the [Runge-Kutta-Fehlberg](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method) (RKF45) method to solve the differential equations, rather than for example 4th order Runge-Kutta or some implicit scheme. This is mostly because RKF45 has a particularly useful way of estimating the numerical error introduced by the integration scheme, and using this estimated error to choose a new appropriate step size. This adaptive step size becomes small when the derivatives of our variables are large, so it very naturally takes care of the coordinate singularities. In particular, we completely get rid of the otherwise present artifact at the vertical axis when using RKF45. 

##### RKF45
 
For some generic first order ordinary differential equation $\frac{dy}{dt} = f(t, y)$, the RKF45 scheme can be written in the following form

$$
\begin{equation}
    k_n = f\left(y + \sum_{i = 1}^{n - 1} B_{ni}k_i\right).
\end{equation}
$$

where the coefficients in the matrix $B$ are given in the article linked to above. The approximation to the function value at time $t + h$, given the function value at $t$ is then 

$$
\begin{equation}
  y(t + h) = y(t) + CH_1 k_1 + CH_2 k_2 + CH_3 k_3 + CH_4 k_4 + CH_5 k_5 + CH_6 k_6,
\end{equation}
$$

where the coefficient in the vector $CH$ are also given in the article. The error estimate associated with the scheme is 

$$
\begin{equation}
  TE = \| CT_1 k_1 + CT_2 k_2 + CT_3 k_3 + CT_4 k_4 + CT_5 k_5 + CT_6 k_6 \|.
\end{equation}
$$

and the new step size is 

$$
\begin{equation}
  h_\text{new} = 0.9 h \left(\frac{\varepsilon}{TE}\right)^{1/5}
\end{equation}
$$

where $\varepsilon$ is the error tolerance we choose. 

### Project structure

- `metric.cpp`: Defines the base `Metric` class from which all other metric classes derive. Defines the methods for implementing the RKF45 integration scheme, coordinate transformations between Cartesian and spherical/Boyer-Lindquist coordinates, function for computing the $\mu = 0$ component of a four-vector given its $\mu = 1, 2, 3$ components. This `Metric` class can also be used directly, in which case we are desccribing Minkowski spacetime. 

- `kerr.cpp`: Defines the `Kerr` class which derives from the `Metric` class. This specializes a lot of the functions from the `Metric` class to Kerr spacetime. This includes the components of the metric tensor, the right-hand side (rhs) of the geodesic equation and the coordinate transformation functions.

- `schwarzschild.cpp`: Completely analogous to the `Kerr` class, the `Schwarzschild` class also derives from the `Metric` class, and describes Schwarzschild spacetime. 

- `scene.cpp`: Defines the `Scene` class, which has methods for initializing and simulating the camera rays. 

- `colorCalculator.cpp`: Defines the `ColorCalculator` class, which has methods for computing the RGB color of a blackbody at a given temperature $T$. 

Instead of defining the position of the camera, its direction, the parameters of the metric etc. inside the code, I use `.json` files to specify all of this. I call these *scene files*, and a given scene fully specifies the state of the simulation. In `src/main.cpp` I read in a scene file and use it to define all of the simulation settings. Although this clutters `main.cpp` it makes it significantly easier to change between simulations, because I can just save different simulation settings in different scene files, and change between them in `main.cpp` by just changing which scene file to read. I have also crated a `scene_template.json` file which shows the structure the code is expecting the scene files to take, together with a description of what each of the parameters describes. 

The `main.cpp` code then initializes a `scene` object from `scene.cpp` (slightly confusing, because this `scene` is not directly related to the scene files themselves. Should change name), which has methods for initializing the rays from the camera and solving the geodesic equation for each of these. `Scene::simulate_camera_rays` contains the main render loop. The generated image is then saved to file in `main.cpp`. 

---

The code is built using `cmake build` and run using `./out/build/Relativistic-renderer/Relativistic-renderer.exe`. 