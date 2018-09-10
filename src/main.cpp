#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

const double Lf = 2.67;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;

  double delta;
  double acceleration;

  h.onMessage([&mpc, &delta, &acceleration](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];//the trajectory we want to follow : x
          vector<double> ptsy = j[1]["ptsy"];//the trajectory we want to follow : y
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];



		 delta = j[1]["steering_angle"];
		 acceleration = j[1]["throttle"];
		 // in order to come over latency, reference: https://discussions.udacity.com/t/how-to-incorporate-latency-into-the-model/257391/3

		  // update state to 100ms in the future for latency, another 10ms for calculation of MPC.
		  double latency = 0.11;
		  px = px + v*cos(psi)*latency;
		  py = py + v*sin(psi)*latency;
		  psi = psi - v*delta / Lf*latency;
		  v = v + acceleration*latency;

		  //transform the coordinate into car's view, reference : Project Q&A: https://youtu.be/bOQuhpz3YfU
		  for (int i = 0; i<int(j[1]["ptsx"].size()); i++) {
			  double gx = (j[1]["ptsx"][i]);
			  double gy = (j[1]["ptsy"][i]);
			  ptsx[i] = (gx - px)*cos(0 - psi) - (gy - py)*sin(0 - psi);
			  ptsy[i] = (gx - px)*sin(0 - psi) + (gy - py)*cos(0 - psi);
		  };
		  
		  //construct a Eigen::VectorXd to get transformed coordinate
		  double* ptrx = &ptsx[0];
		  Eigen::Map<Eigen::VectorXd> ptsx_transform(ptrx, 6);

		  double* ptry = &ptsy[0];
		  Eigen::Map<Eigen::VectorXd> ptsy_transform(ptry, 6);

          /*
          * TODO: Calculate steering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */

		  //========================= step1.fit the trajectory we want to follow ===========================//
		  Eigen::VectorXd coeffs = polyfit(ptsx_transform, ptsy_transform, 3); // free to adjust order
		  //========================= step1.fit the trajectory we want to follow ===========================//

		  //====================== step2.calculate cross track error and psi error =========================//
		  double cte = polyeval(coeffs, 0);
		  double epsi = psi - atan((3 * coeffs[3] * px * px + 2 * coeffs[2] * px + coeffs[1])); // depend on order : derivative
		  
		  //====================== step2.calculate cross track error and psi error =========================//

		  //================================== step3.apply MPC method ======================================//
		  //  needn't to use iteration? no, we don't use iteration
		  
		  //Display the MPC predicted trajectory 

		  vector<double> mpc_x_vals;
		  vector<double> mpc_y_vals;

		  //Display the waypoints/reference line
		  vector<double> next_x_vals;
		  vector<double> next_y_vals;
		  
		  Eigen::VectorXd state(6);
		  

		  state << 0, 0, 0, v, cte, epsi;

		  auto vars = mpc.Solve(state, coeffs); //we need to get "a" and "delta" after many steps
		  
		  //reference : Project Q&A: https://youtu.be/bOQuhpz3YfU
		  double poly_distance = 2.5;
		  int num_points = 25;
		  for (int i = 1; i < num_points; i++)
		  {
			  next_x_vals.push_back(poly_distance*i);
			  next_y_vals.push_back(polyeval(coeffs, poly_distance*i));
		  }
		  
		  for (int i = 2; i < int(vars.size()); i++)
		  {
			  if (i % 2 == 0)
			  {
				 mpc_x_vals.push_back(vars[i]);
			  }
			  else
			  {
				  mpc_y_vals.push_back(vars[i]);
			  }
		  }
		  
		  //================================== step3.apply MPC method ======================================//
          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
		  const double Lf = 2.67;
		  double steer_value = vars[0]/(deg2rad(25) * Lf);
		  double throttle_value = vars[1];
		  
		  msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

		  //delta = steer_value;
		  //acceleration = throttle_value;

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;



          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
