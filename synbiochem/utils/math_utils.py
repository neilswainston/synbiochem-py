'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-member
import gurobipy


def isclose(value_1, value_2, rel_tol=1e-09, abs_tol=0.0):
    '''Compares floating point numbers.'''
    return abs(value_1 - value_2) <= max(rel_tol * max(abs(value_1),
                                                       abs(value_2)), abs_tol)


def linprog(c_vector, a_eq, b_eq, bounds):
    '''Solve linear programming problem.'''
    model = gurobipy.Model()
    model.Params.OutputFlag = 0  # Silent mode

    # Create variables:
    for obj, bound in zip(c_vector, bounds):
        model.addVar(bound[0], bound[1], obj, gurobipy.GRB.CONTINUOUS)

    # Integrate new variables:
    model.update()

    # Add constraints:
    for aaa, bbb in zip(a_eq, b_eq):
        model.addConstr(gurobipy.LinExpr(aaa, model.getVars()),
                        gurobipy.GRB.EQUAL,
                        bbb)

    model.optimize()

    return model.status == gurobipy.GRB.status.OPTIMAL, \
        [var.X for var in model.getVars()] \
        if model.status == gurobipy.GRB.status.OPTIMAL else None
