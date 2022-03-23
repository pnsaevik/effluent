from effluent import modeldata


class Test_from_profile:
    def test_returns_dataset_with_standardized_names(self):
        m = modeldata.Modeldata.from_single_profile(
            z=[0, 10, 50],
            u=[.3, .2, .1],
            v=[.03, .02, .01],
            dens=[1023, 1024, 1025],
        )

        p = m.profile()

        assert set(p.variables) == {'z', 'dens', 'u', 'v', 'w', 'lat', 'lon', 'date'}
