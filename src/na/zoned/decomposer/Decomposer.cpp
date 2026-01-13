//
// Created by cpsch on 11.12.2025.
//
#include "na/zoned/decomposer/decomposer.hpp"

//#include "ir/operations/StandardOperation.hpp"
#include "ir/operations/CompoundOperation.hpp"

#include <complex>
#include <vector>

namespace na::zoned {
Decomposer::Decomposer(int n_qubits) {
  N_qubits=n_qubits;
};


auto Decomposer::combine_quaternions(const std::array<qc::fp,4>& q1,
                         const std::array<qc::fp,4>& q2) -> std::array<qc::fp,4> {
  std::array<qc::fp,4> new_quat{};
  new_quat[0]=q1[0]*q2[0]-q1[1]*q2[1]-q1[2]*q2[2]-q1[3]*q2[3];
  new_quat[1]=q1[0]*q2[1]+q1[1]*q2[0]+q1[2]*q2[3]-q1[3]*q2[2];
  new_quat[2]=q1[0]*q2[2]-q1[1]*q2[3]+q1[2]*q2[0]+q1[3]*q2[1];
  new_quat[3]=q1[0]*q2[3]+q1[1]*q2[2]-q1[2]*q2[1]+q1[3]*q2[0];
  return new_quat;
}

auto Decomposer::convert_gate_to_quaternion(std::reference_wrapper<const qc::Operation> op) -> std::array<qc::fp,4> {
  assert(op.get().getNqubits() == 1);
  std::array<qc::fp,4> quat;
  //TODO: Figure out if I need -i factor
  //TODO: Phase?
  if (op.get().getType() == qc::RZ || op.get().getType() == qc::P) {
        quat={cos(op.get().getParameter().front()),0,0,sin(op.get().getParameter().front())};
      } else if (op.get().getType() == qc::Z) {
        quat={0,0,0,1};
      } else if (op.get().getType() == qc::S) {
        quat={cos(qc::PI_4),0,0,sin(qc::PI_4)};
      } else if (op.get().getType() == qc::Sdg) {
        quat={cos(-qc::PI_4),0,0,sin(-qc::PI_4)};
      } else if (op.get().getType() == qc::T) {
        quat={cos(qc::PI_4/2),0,0,sin(qc::PI_4/2)};
      } else if (op.get().getType() == qc::Tdg) {
        quat={cos(-qc::PI_4/2),0,0,sin(-qc::PI_4/2)};
      } else if (op.get().getType() == qc::U) {
        quat=combine_quaternions(combine_quaternions({cos(op.get().getParameter().at(1)/2),0,0,sin(op.get().getParameter().at(1)/2)}
          ,{cos(op.get().getParameter().front()/2),0,sin(op.get().getParameter().front()/2),0}),
            {cos(op.get().getParameter().at(2)/2),0,0,sin(op.get().getParameter().at(2)/2)});
        } else if (op.get().getType() == qc::U2) {
          quat=combine_quaternions(combine_quaternions({cos(op.get().getParameter().front()/2),0,0,sin(op.get().getParameter().front()/2)}
          ,{cos(qc::PI_2),0,sin(qc::PI_2),0}),
            {cos(op.get().getParameter().at(1)/2),0,0,sin(op.get().getParameter().at(1)/2)});
        } else if (op.get().getType() == qc::RX) {
          quat={cos(op.get().getParameter().front()/2),sin(op.get().getParameter().front()/2),0,0};
        } else if (op.get().getType() == qc::RY) {
          quat={cos(op.get().getParameter().front()/2),0,sin(op.get().getParameter().front()/2),0};
        } else if (op.get().getType() == qc::H) {
          //TODO: Double check this
          quat=combine_quaternions(combine_quaternions({1,0,0,0}
          ,{cos(qc::PI_4),0,sin(qc::PI_4),0}),
            {cos(qc::PI_2),0,0,sin(qc::PI_2)});
        } else if (op.get().getType() == qc::X) {
          quat={0,1,0,0};
        } else if (op.get().getType() == qc::Y) {
          quat={0,1,0,0};
        } else if (op.get().getType() == qc::Vdg) {
          quat=combine_quaternions(combine_quaternions({cos(qc::PI_4),0,0,sin(qc::PI_4)}
          ,{cos(-qc::PI_4),0,sin(-qc::PI_4),0}),
            {cos(-qc::PI_4),0,0,sin(-qc::PI_4)});
        } else if (op.get().getType() == qc::SX) {
          quat=combine_quaternions(combine_quaternions({cos(-qc::PI_4),0,0,sin(-qc::PI_4)}
            ,{cos(qc::PI_4),0,sin(qc::PI_4),0}),
          {cos(qc::PI_4),0,0,sin(qc::PI_4)});
        } else if (op.get().getType() == qc::SXdg ||
                   op.get().getType() == qc::V) {
          quat=combine_quaternions(combine_quaternions({cos(-qc::PI_4),0,0,sin(-qc::PI_4)}
            ,{cos(-qc::PI_4),0,sin(-qc::PI_4),0}),
            {cos(qc::PI_4),0,0,sin(qc::PI_4)});
        } else {
          // if the gate type is not recognized, an error is printed and the
          // gate is not included in the output.
          std::ostringstream oss;
          oss << "\033[1;31m[ERROR]\033[0m Unsupported single-qubit gate: "
              << op.get().getType() << "\n";
          throw std::invalid_argument(oss.str());
        }
        return quat;
      }

  auto Decomposer::get_U3_angles_from_quaternion(std::array<qc::fp,4> quat) -> std::array<qc::fp,3> {
  //TODO: Is there a prescribed eps somewhere? Else where should THIS eps be defined
  qc::fp eps=1e-5;
  qc::fp theta;
  qc::fp phi;
  qc::fp lambda;
  if (abs(quat[0])>eps||abs(quat[3])>eps) {
    theta=2.*std::atan2(std::sqrt(quat[2]*quat[2]+quat[1]*quat[1]),std::sqrt(quat[0]*quat[0]+quat[3]*quat[3]));
    qc::fp alph_1=std::atan2(quat[3],quat[0]); //phi+ lambda
    if (abs(quat[1])>eps||abs(quat[2])>eps) {
      qc::fp alph_2=-1*std::atan2(quat[1],quat[2]);
      phi=alph_1+alph_2; //phi
      lambda=alph_1-alph_2;
    } else {
      //TODO: Which was set to zero???
      phi=0;
      lambda=2*alph_1;
    }
  } else {
    theta=qc::PI; //or § PI if sin(theta/2)=-1... Relevant? Or is the -1= exp(i*pi) just global phase
    if (abs(quat[1])>eps||abs(quat[2])>eps) {
      phi=0;
      lambda= 2*std::atan2(quat[1],quat[2]);
      //2*tan(q1/q2) = phi-lambda, but can't be disentangled
    } else {
      //This should never happen! Exception??
      phi=0.;
      lambda=0.;
    }
  }
  return {theta,phi,lambda};
}


 auto Decomposer::calc_theta_max(const SingleQubitGateLayer& layer) -> qc::fp {
  //TODO: Error Handling improvement?
  qc::fp theta_max=0;
  for (auto  gate :layer) {
    if (gate.get().getType()== qc::U){
    std::vector<qc::fp> params=gate.get().getParameter();
    if (abs(params[0])>theta_max) {theta_max=abs(params[0]);}
    } else {
    // if the gate type is not U3, an error is printed.
    std::ostringstream oss;
    oss << "\033[1;31m[ERROR]\033[0m Unsupported single-qubit gate: "
        << gate.get().getType() << "\n";
    throw std::invalid_argument(oss.str());
    }
  }
    return theta_max;
}

//TOD0:Return vector of structs
auto Decomposer::transform_to_U3(const std::vector<SingleQubitGateLayer>& layers) const
    -> std::vector<SingleQubitGateLayer> {
  //How to sort gates with same qubit
  //How to get Matrix Representations??
  std::vector<SingleQubitGateLayer> new_layers;
  for (const auto& layer: layers) {
    std::vector<std::vector<std::reference_wrapper<const qc::Operation>>> gates(this->N_qubits);
    SingleQubitGateLayer new_layer;
    for (auto gate:layer) {
      gates[gate.get().getNtargets()].push_back(gate);
    }

    for (auto qubit_gates:gates) {
      if (!qubit_gates.empty()) {

        std::array<qc::fp,4> quat = convert_gate_to_quaternion(qubit_gates[0]);
        for (auto i=1;i<qubit_gates.size();i++) {
          quat=combine_quaternions(quat,convert_gate_to_quaternion(qubit_gates[i]));
        }
        std::array<qc::fp,3> angles=get_U3_angles_from_quaternion(quat);
        //TODO: Don't return Op here?
        const qc::Operation *new_op=new qc::StandardOperation(qubit_gates[0].get().getTargets().front(),qc::U,{angles[0],angles[1],angles[2]});
        new_layer.emplace_back(*new_op);
      }
    }
    new_layers.push_back(new_layer);
  }
  return new_layers;
}


auto Decomposer::decompose(
    const std::pair<std::vector<SingleQubitGateLayer>,
                    std::vector<TwoQubitGateLayer>>& schedule)
    -> std::pair<std::vector<SingleQubitGateLayer>,
                 std::vector<TwoQubitGateLayer>> {
  //TODO:make U3Layer vector<struct ...>
  std::vector<SingleQubitGateLayer> U3Layers=this->transform_to_U3(schedule.first);
  std::vector<SingleQubitGateLayer> NewSingleQubitLayers=std::vector<SingleQubitGateLayer>{};

  for (const auto& layer: U3Layers) {
    qc::fp theta_max=calc_theta_max(layer);
    SingleQubitGateLayer FrontLayer;
    SingleQubitGateLayer MidLayer;
    SingleQubitGateLayer BackLayer;
    SingleQubitGateLayer NewLayer;


    for (auto gate:layer){
      //TODO: ONly ANgles at this point
      std::vector<qc::fp> params=gate.get().getParameter();
      qc::fp theta=params[0];
      qc::fp phi_minus=params[1];
      qc::fp phi_plus=params[2];
      //TODO: Decide whether to make this a separate function for each gate in rel to theta_max
      //TODO: Add mod 2PI/TAU
      qc::fp alpha, chi, beta;
      if (theta==theta_max) {
        alpha=qc::PI_2;
        chi=qc::PI;
      }else {
        qc::fp kappa=sqrt((sin(theta/2)*sin(theta/2))/(sin(theta_max)*sin(theta_max)-(sin(theta/2)*sin(theta/2))));
        alpha=atan(cos(theta_max/2)*kappa);
        chi=fmod(2*atan(kappa),qc::TAU);
      }
      if (theta!=0) {
        beta=theta_max<0?qc::PI_2:-1*qc::PI_2;
      }else {
         beta=0;
      }
      qc::fp gamma_plus=fmod(phi_plus-(alpha+beta),qc::TAU);
      qc::fp gamma_minus=fmod(phi_minus-(alpha-beta),qc::TAU);

      //U3(theta,phi_min(phi),phi_plus(lamda)->Rz(gamma_minus)GR(theta_max/2, PI_2)Rz(chi)GR(-theta_max/2,PI_2)RZ(gamma_plus)
      // GR(theta_max/2, PI_2)==Global Y due to PI_2


      //TODO: Does this work??
      auto sop=qc::StandardOperation(gate.get().getTargets().front(),qc::RZ,{gamma_minus});
      qc::Operation *op=&sop;
      FrontLayer.emplace_back(std::reference_wrapper<qc::Operation>(*op));

      sop=qc::StandardOperation(gate.get().getTargets().front(),qc::RZ,{chi});
      op=&sop;
      MidLayer.emplace_back(std::reference_wrapper<qc::Operation>(*op));

      sop=qc::StandardOperation(gate.get().getTargets().front(),qc::RZ,{gamma_plus});
      op=&sop;
      MidLayer.emplace_back(std::reference_wrapper<qc::Operation>(*op));
    }//gate::layer



    std::vector< std::unique_ptr<qc::Operation>> GR_plus;
    std::vector< std::unique_ptr<qc::Operation>>  GR_minus;

    for (auto i=0; i<this->N_qubits; ++i) {
      GR_plus.emplace_back(new qc::StandardOperation(i,qc::RY,{theta_max/2}));
      GR_minus.emplace_back(new qc::StandardOperation(i,qc::RY,{-1*theta_max/2}));
    }

    for (auto gate:FrontLayer) {
      NewLayer.push_back(gate);
    }

    auto cop= qc::CompoundOperation(std::move(GR_plus),true);
    qc::Operation *ryp=&cop;
    NewLayer.emplace_back(*ryp);

    for (auto gate:MidLayer) {
      NewLayer.push_back(gate);
    }
    cop=qc::CompoundOperation(std::move(GR_minus),true);
    qc::Operation *rym=&cop;
    NewLayer.emplace_back(*rym);

    for (auto gate:BackLayer) {
      NewLayer.push_back(gate);
    }
    NewSingleQubitLayers.push_back(NewLayer);
  }//layer::SingleQubitLayers
  return std::pair{NewSingleQubitLayers,schedule.second};
}
}// namespace na::zoned