
    """
    pca = PCA(n_components=2)
    pca.fit(coordinates_rotated[1])
    X_pca = pca.transform(coordinates_rotated[1])
    X_new = pca.inverse_transform(X_pca)

    print("original shape:   ", coordinates_rotated[1].shape)
    print("transformed shape:", X_pca.shape)

    pca.fit(coordinates_rotated[0])
    Y_pca = pca.transform(coordinates_rotated[0])
    Y_new = pca.inverse_transform(Y_pca)

    pca.fit(coordinates_rotated[2])
    Z_pca = pca.transform(coordinates_rotated[2])
    Z_new = pca.inverse_transform(Z_pca)


    #plt.scatter(coordinates_rotated[1][:, 0], coordinates_rotated[1][:, 1], alpha=0.2)
    plt.scatter(X_new[:, 0], X_new[:, 1], alpha=0.8)
    plt.scatter(Y_new[:, 0], Y_new[:, 1], alpha=0.2)
    plt.scatter(Z_new[:, 0], Z_new[:, 1], alpha=0.5)

    plt.axis('equal')
    plt.show()
    """
    """

        pca = PCA(n_components=2)
        for molecule in coordinates_rotated:
            pca.fit(molecule)
            transformed_pca = pca.transform(molecule)
            inverse_transformed = pca.inverse_transform(transformed_pca)
            plt.scatter(inverse_transformed[:, 0], inverse_transformed[:, 1], alpha=0.8)
            plt.axis('equal')
        plt.show()
    """
