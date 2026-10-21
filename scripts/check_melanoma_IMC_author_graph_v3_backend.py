import ast,importlib.util,inspect,json,textwrap,unittest
from pathlib import Path
root=Path(__file__).resolve().parents[1]
def load(path,name):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
fixtures=load(root/'tests/test_melanoma_IMC_complete_author_graph_states.py','author_v2_fixtures')
def setup(cls):
 cls.m=load(root/'scripts/analyze_melanoma_IMC_author_graph_states_v3.py','actual_author_core_v2')
 ordinary=cls.m.analyze
 def fixture(bundle,parameters):
  bundle['author_neighbor_edges']={roi:bundle['native_knn_edges'](bundle['native_coordinates'][idx],10) for roi,idx in bundle['frame'].groupby('sample_id').indices.items()}
  def forbidden(*args):raise AssertionError('native kNN not permitted')
  bundle['native_knn_edges']=forbidden
  return ordinary(bundle,parameters)
 cls.m.analyze=fixture
fixtures.AuthorAdjacencyMath.setUpClass=classmethod(setup)
method=fixtures.AuthorAdjacencyMath.test_full_synthetic_wiring_2000_draws_and_18_opportunities
source=textwrap.dedent(inspect.getsource(method)).replace('private_frozen_step_author_node_native_geometry_v1_source_summary.json','private_author_node_counter_qualification.json').replace('COMPLETE_AUTHOR_NODE_BIJECTION_AND_EDGE_PREFIX_NATIVE_COORDINATE_REVIEW','COMPLETE_EXPLICIT_COUNTER_BASE_AND_LABEL_BIJECTION')
namespace=method.__globals__.copy();exec(compile(source,'<updated proof-path fixture>','exec'),namespace);fixtures.AuthorAdjacencyMath.test_full_synthetic_wiring_2000_draws_and_18_opportunities=namespace[method.__name__]
def streamsetup(cls):cls.m=load(root/'scripts/measure_melanoma_IMC_complete_author_graph_states_v3.py','actual_author_stream_v2')
fixtures.CompleteAuthorEdges.setUpClass=classmethod(streamsetup)
keys=load(root/'tests/test_melanoma_IMC_explicit_counter_bridge.py','actual_reversible_key_contracts')
v1=(root/'scripts/analyze_melanoma_IMC_author_graph_states_v1.py').read_text();v2=(root/'scripts/analyze_melanoma_IMC_author_graph_states_v3.py').read_text();assert v1==v2.replace('private_author_node_counter_qualification.json','private_frozen_step_author_node_native_geometry_v1_source_summary.json').replace('COMPLETE_EXPLICIT_COUNTER_BASE_AND_LABEL_BIJECTION','COMPLETE_AUTHOR_NODE_BIJECTION_AND_EDGE_PREFIX_NATIVE_COORDINATE_REVIEW')
def functions(path):return {x.name:ast.dump(x,include_attributes=False) for x in ast.parse(path.read_text()).body if isinstance(x,(ast.FunctionDef,ast.ClassDef))}
a=functions(root/'scripts/measure_melanoma_IMC_complete_author_graph_states_v1.py');b=functions(root/'scripts/measure_melanoma_IMC_complete_author_graph_states_v3.py');assert all(a[name]==b[name] for name in ['DigestReader','resolve_rows','selected_edges'])
loader=unittest.TestLoader();suite=unittest.TestSuite([loader.loadTestsFromModule(fixtures),loader.loadTestsFromModule(keys)]);result=unittest.TextTestRunner(verbosity=2).run(suite);assert result.wasSuccessful() and result.testsRun==22 and not result.skipped
print(json.dumps(dict(status='CURRENT_AUTHOR_GRAPH_REVERSIBLE_BRIDGE_TESTS_PASSED',tests=22,skipped=0,all_statistical_math_unchanged=True,three_graph_stream_functions_AST_identical=True)))
