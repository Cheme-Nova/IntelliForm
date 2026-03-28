import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.memory import memory

class IntelliFormController:
    def __init__(self):
        from modules.optimizer import run_optimization
        from modules.pareto import run_pareto_optimization
        from modules.bayesian import run_bayesian_optimization
        from modules.eco import compute_ecometrics
        from modules.qsar import predict_properties
        from modules.regulatory import get_blend_report, generate_vertical_regulatory_report
        from modules.stability import predict_stability
        from modules.carbon import calculate_carbon_credits
        from modules.certifications import run_certification_oracle
        from modules.agents import run_agent_swarm
        from modules.nlp import parse_request

        self.parse_request = parse_request
        self.run_optimization = run_optimization
        self.run_pareto = run_pareto_optimization
        self.run_bayesian = run_bayesian_optimization
        self.compute_eco = compute_ecometrics
        self.predict_qsar = predict_properties
        self.get_reg = get_blend_report
        self.get_vreg = generate_vertical_regulatory_report
        self.predict_stability = predict_stability
        self.calc_carbon = calculate_carbon_credits
        self.run_cert = run_certification_oracle
        self.run_agents = run_agent_swarm
        self.memory = memory

    def run(self, input_text, vertical, batch_size, opt_mode, constraints, db):
        parsed = self.parse_request(input_text)
        result = self.run_optimization(db, constraints, opt_mode)
        eco = self.compute_eco(result.blend, db)
        reg = self.get_reg(result.blend)
        vreg = self.get_vreg(result.blend, db, vertical)
        stability = self.predict_stability(result.blend, db)
        carbon = self.calc_carbon(result.blend, db, batch_size)
        cert = self.run_cert(result.blend, db, vertical, result.bio_pct)
        agents = self.run_agents(result, parsed)
        self.memory.record("formulation_generated", result.blend, vertical)
        return {
            "parsed": parsed,
            "result": result,
            "eco": eco,
            "reg": reg,
            "vreg": vreg,
            "stability": stability,
            "carbon": carbon,
            "cert": cert,
            "agents": agents
        }

controller = IntelliFormController()
