def test_webapp_imports():
    import webapp
    assert hasattr(webapp, "__version__")
