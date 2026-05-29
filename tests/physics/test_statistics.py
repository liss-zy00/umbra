import numpy as np
from src.physics.diagnostics.statistics import edge_radial_mean, pearson_correlation_coefficient

def test_edge_radial_mean():
    # Array de 10 elementos. Queremos os 20% finais (2 elementos),
    # mas o 'min_points' padrão é 3. Portanto, deve pegar os últimos 3.
    values = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    
    result = edge_radial_mean(values, fraction=0.2, min_points=3)
    
    # Média de [8.0, 9.0, 10.0] = 9.0
    assert np.isclose(result, 9.0)

def test_edge_radial_mean_with_nans():
    values = [1.0, 2.0, 3.0, np.nan, 8.0, 9.0, 10.0, np.nan]
    result = edge_radial_mean(values, fraction=0.5, min_points=2)
    # Valores válidos: [1, 2, 3, 8, 9, 10]. Últimos 3 (50% de 6): [8, 9, 10]
    assert np.isclose(result, 9.0)

def test_pearson_correlation_coefficient():
    a = np.array([1.0, 2.0, 3.0, 4.0])
    b = np.array([2.0, 4.0, 6.0, 8.0]) # Correlação perfeita (1.0)
    
    result = pearson_correlation_coefficient(a, b)
    assert np.isclose(result, 1.0)

def test_pearson_invalid_data():
    a = np.array([1.0, 1.0, 1.0, 1.0]) # Desvio padrão zero
    b = np.array([2.0, 4.0, 6.0, 8.0])
    assert pearson_correlation_coefficient(a, b) is None